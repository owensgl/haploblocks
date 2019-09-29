#!/bin/bash

function usage ()
{
    echo "Usage: $(basename "$0") [--tmpdir TMPDIR] [--max-mem-mb MEMMB] PARAMSFILE OUTSNP OUTINDEL INPUTVCF+

    Selects good variants from the raw input vcfs provided, and writes
    the output to OUTSNP (snps) and OUTINDEL (indels).

      If INPUTVCF starts with '@', then the parameter is assumed to a
      filename containing a list of vcf input filenames.

    PARAMSFILE is a file describing the filtering job to apply. it is
    a linebase file with a series of directives. Each line can be:

      #SM <samplename>          # listing of the sample names covered
                                # in the input. one such line per sample.

      #MAX_AN <num>             # the number of chromosomes covered in the
                                # input. typically 2*num_samples (diploid).

      #FILTER <expr>            # a filtering expression, using $V[stat] to
                                # access VCF statistics columns. Those are
                                # embedded into an AWK script.
"

}

set -eo pipefail

# shut up perl.  the LANG environment variable is propagated from the
# login node into the compute job. but a different set of locales are
# installed.
export LANG=C

SCRIPTSDIR="$(readlink -f $(dirname "$0"))"
VCF2STATS_BIN="${SCRIPTSDIR}/vcf2stats.pl"
VCFFILTER_BIN="${SCRIPTSDIR}/../filters/test.awk)"
SELECTSITES_BIN="${SCRIPTSDIR}/vcf2selectsites.pl"
RESERVOIR_BIN="${SCRIPTSDIR}/reservoir.py"
BGZIP_BIN="${CONDA_PREFIX}/bin/bgzip"
PARALLEL_BIN="${CONDA_PREFIX}/bin/parallel"
AWK_BIN="${CONDA_PREFIX}/bin/awk"

#
# Number of samples to randomly select, when
# sampling the distribution
#
RESERVOIR_SAMPLES=1000000

CFGTMP=/tmp
MAX_MEM_MB=$((4*1024))
INPUTS=()
OUTPUTS=()
PARAM_FILE=""
REFERENCE=""

while [[ $# -gt 0 ]]; do
    arg="$1"
    shift
    case "$arg" in
	--tmpdir)
	    CFGTMP="$1"
	    shift
	    ;;
	--max-mem-mb)
	    MAX_MEM_MB="$1"
	    shift
	    ;;
	--reference)
	    REFERENCE="$1"
	    shift
	    ;;
	--help|-h)
	    usage
	    exit 0
	    ;;
	@*)
	    # read filenames from file
	    infile="${arg##@}"
	    while read inputvcf rest; do
		INPUTS+=( "$inputvcf" )
	    done < "$infile"
	    ;;
	*)
	    if [[ -z "$PARAM_FILE" ]]; then
		PARAM_FILE="$arg"
	    elif [[ ${#OUTPUTS[@]} -lt 2 ]]; then
		OUTPUTS+=( "$arg" )
	    else
		INPUTS+=( "$arg" )
	    fi
	    ;;
    esac
done

set -x
START_TS=$(printf "%(%s)T" -1)

# Keep some extra for stack
MAX_MEM_MB=$((MAX_MEM_MB - 512))
if [[ ${MAX_MEM_MB} -lt 1024 ]]; then
    echo "Provide more ram to --max-mem-mb <AMOUNT_MB>" >&2
    exit 1
fi
JAVA_OPTIONS="-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx${MAX_MEM_MB}m"

if [[ ${#OUTPUTS[@]} -lt 2 ]]; then
    echo "missing output(s)" >&2
    usage >&2
    exit 1
fi

if [[ ${#INPUTS[@]} -lt 1 ]]; then
    echo "Missing input file(s)" >&2
    usage >&2
    exit 1
fi

if [[ -n "$REFERENCE" && ! -e "$REFERENCE" ]]; then
    echo "Reference file $REFERENCE is invalid." >&2
    exit 1
fi

XTMP="${SLURM_TMPDIR:-$CFGTMP}" # use local disk if on cluster
mkdir -p "$XTMP" "$XTMP"/gatk
workdir=$(mktemp -d -p $XTMP tmp.midas.XXXXXX);
mkdir -p "${workdir}/gatk"
TMPFILES=()
trap 'rm --one-file-system -rf -- ${workdir}; if [[ ${#TMPFILES[@]} -gt 0 ]]; then rm -f --one-file-system -- "${TMPFILES[@]}"; fi' EXIT;

# Prefixes each line of stdin with a name and
# the time delta in seconds since the first line was read
function prefix_log () { # PREFIX [display_time]
    (
	set +x;
	local start_ts=$(printf "%(%s)T" -1)
	local prefix="$1"
	local do_time="${2:-""}"
	local LINE
	if [[ -n "${do_time}" ]]; then
	    while IFS="\n" read LINE; do
		printf "[%s %5d] %s\n" "$prefix" "$(( $(printf "%(%s)T" -1) - start_ts))" "$LINE"
	    done
	    printf "[%s %5d] %s\n" "$prefix" "$(( $(printf "%(%s)T" -1) - start_ts))" "<EOF>"
	else
	    while IFS="\n" read LINE; do
		printf "[%s] %s\n" "$prefix" "$LINE"
	    done
	    printf "[%s] %s\n" "$prefix" "<EOF>"
	fi
    )
}


function calculate_distribution ()
{ (
    set -eo pipefail
    local statname="$1";

    {
	# Only compute statistics on rows with FILTER==PASS.
	# Those rows are the only ones considered by downstream
	# filters.
	"${AWK_BIN}" -vOFMT=%20.7f \
'BEGIN { i=-1; filtercol=-1; }
NR==1 {
    filtercol=-1;
    for (i=1; i<=NF; i++) {
        if (tolower($i) == "filter") { filtercol = i; }
    }
    if (filtercol == -1) {
        print "Error: no FILTER column." > "/dev/stderr"; is_error = 1;
        exit 1;
    }
    print $0; # keep header
}
NR > 1 && NF > 0 && filtercol > 0 {
   if ($filtercol == "PASS") { print $0; }
}
'
    } | {
	# Select N random records.  Avoids numerical errors for large
	# cohorts. Loss of precision with 53 bit floating point happens after
	# millions of entries. (got billions).
	"${RESERVOIR_BIN}" "${RESERVOIR_SAMPLES}" --col "${statname}"
    } | {
	# compute sample mean and stddev. we use the population equation,
	# and then remove bias with Bessel's correction to obtain the
	# sample stddev.
	"${AWK_BIN}" -vOFMT=%20.7f -vstatname="${statname}" '
BEGIN { sum=0; sum_squares=0; N=0; datum=0; statcol=-1; maxval="-inf"+0; minval="+inf"+0; }
NR==1 {
  for (i=1; i<=NF; i++) {
      if (tolower(statname) == tolower($i)) { statcol=i; break; }
  }
  if (statcol == -1) {
     print "Error: no " statname " column." > "/dev/stderr"; is_error = 1; exit 1;
  }
}
NR > 1 && NF {
  N+=1;
  datum=($statcol+0);
  sum+=datum;
  if (datum > maxval) { maxval = datum; }
  if (datum < minval) { minval = datum; }
  sum_squares+=(datum)*(datum); # **2 not portable.
  #print $statcol, sum, sum_squares;
}
END {
  if (is_error == 1) exit 1;
  if (N == 0) { exit 1; }
  mean=sum/N;
  pop_variance=(sum_squares/N)-(mean*mean); # **2 not portable
  pop_stddev=sqrt(pop_variance);
  if (N > 1) {
      sam_variance=(N/(N-1)) * pop_variance; # bessel correction
      sam_stddev=sqrt(sam_variance);
  } else {
      sam_variance = "+nan"+0;
      sam_stddev= "+nan"+0;
  }
  print "N", "sum", "sum_sq", "mean", "pop_var", "pop_sd", "sam_var", "sam_sd", "min", "max";
  print N, sum, sum_squares, mean, pop_variance, pop_stddev, sam_variance, sam_stddev, minval, maxval;
}'
    }
) }

#
# Generate the filtration script.  The same filter script is used for
# all input VCF files.
#
function generate_filter_prog ()
{
	local directive="" rest="" sample_count=0
	local filters=() MAX_AN=""
	local params_file="$1"
	while read directive rest; do
	    case "$directive" in
		"#SM")
		    sample_count=$((sample_count+1))
		    ;;
		"#MAX_AN")
		    MAX_AN="$rest"
		    ;;
		"#FILTER")
		    filters+=("$rest")
		    ;;
		*)
		    echo "WARN: Unknown directive in params file: $directive $rest" >&2
		    ;;
	    esac
	done < "$params_file"
	if [[ -z "$MAX_AN" ]]; then
	    echo "MAX_AN missing from params files ${params_file}" >&2
	    return 1
	fi
	if [[ "${#filters[@]}" -lt 1 ]]; then
	    echo "No filters in ${params_file}" >&2
	    return 1
	fi

	{
	    cat <<EOF
function isnan(x) { return ((x+0) == "+nan"+0); } # awk has its own rules
function is_pass(V, MAX_AN, MEAN_DP, SD_DP)
{
        # BEGIN USER SECTION
$(for ((i=0; i<${#filters[@]}; i++)); do
      filter="${filters[i]}";
      printf 'if (!(%s)) { filter_step[%d]++; return 0;}\n' "$filter" "$i";
done)
        # END   USER SECTION
        return 1; # passed
}
BEGIN { MAX_AN=${MAX_AN}+0;
        error_exit = 0;
        if (MEAN_DP=="" || SD_DP=="") { print "Error: variables MEAN_DP and SD_DP must be set." > "/dev/stderr"; error_exit = 1; exit 1; }
        MEAN_DP=MEAN_DP+0;
        SD_DP=SD_DP+0;
        pass_count = 0;
        printf("Filter Globals: MAX_AN=%9.6f MEAN_DP=%9.6f SD_DP=%9.6f\n", MAX_AN, MEAN_DP, SD_DP) > "/dev/stderr";
$(for ((i=0; i<${#filters[@]}; i++)); do
      filter="${filters[i]}";
      printf '        filter_step[%d] = 0;\n' "$i";
done)

      }
NR==1 {
    split("Chr|Pos|Type|Alt_Length|Filter|Qual|Het_Percent|AN|BaseQRankSum|ClippingRankSum|DP|ExcessHet|FS|InbreedingCoeff|MQ|MQRankSum|QD|ReadPosRankSum|SOR|AF_0|AF_1|AF_2", col_names, "|")
    for (colid in col_names) {
	colname = col_names[colid];
	for (i=1; i<=NF; i++) {
	    if (tolower(colname) == tolower(\$i)) { V[colname] = i; }
	}
	if (!V[colname]) { print "Error: Missing column: " colname > "/dev/stderr"; exit 1; }
    }
    print "Chr\tPos\tType"; #header
}
NR > 1 { if (is_pass(V, MAX_AN, MEAN_DP, SD_DP)) {
        pass_count++;
        print \$V["Chr"] "\t" \$V["Pos"] "\t" \$V["Type"]; }
}

END {
    if (error_exit) { exit error_exit; }
    print  "Filter results:" > "/dev/stderr";
    total_input = NR - 1;
    printf("  total input   : %15d\n", total_input) > "/dev/stderr";
$(for ((i=0; i<${#filters[@]}; i++)); do
    printf "    filter_pct = 0.0;\n"
    printf "    if (total_input > 0) { filter_pct = filter_step[%d] * 100.0 / total_input; }\n" "$i"
    printf "    printf(\"  failed step %2d: %%15d (%%2.1f%%%%)\\\n\", filter_step[%d], filter_pct) > \"/dev/stderr\";\n" "$i" "$i";
done)
    filter_pct = 0.0;
    if (total_input > 0) { filter_pct = pass_count * 100.0 / total_input; }
    printf("  total pass    : %15d (%2.1f%%)\n", pass_count, filter_pct) > "/dev/stderr";
}
EOF
	}
}

STATS_FILES=()
PASS_FILES=()
for ((i=0; i<${#INPUTS[@]}; i++)); do
    STATS_FILES[i]=$(printf "%s.%05d.stats" "${workdir}/vcf" $i);
    PASS_FILES[i]=$(printf "%s.%05d.pass" "${workdir}/vcf" $i);
    GOLD_SNPS[i]=$(printf "%s.%05d.snps.vcf.gz" "${workdir}/gold" $i);
    GOLD_INDELS[i]=$(printf "%s.%05d.indels.vcf.gz" "${workdir}/gold" $i);
done

# Extract just the columns needed. Run in parallel.
function vcf_input_stream ()
{
    (
	set +x;
	local invcf
	local count=0 start=0 end=0
	for invcf in "$@"; do
	    count=$((count+1))
	    start=$(date +%s)
	    echo "[$start] Feeding VCF $count/$#: $invcf" >&2
	    if [[ ! -e "$invcf" ]]; then
		echo "invalid vcf file $invcf" >&2
		return 1
	    fi
	    zcat "$invcf" || {
		echo "decompression of $invcf failed (exit $?)" >&2
		return 1
	    }
	    echo "Took $(($(date +%s) - start)) seconds." >&2
	done
    )
}

(

    :
    : Extract stats colums
    :
    function vcf_stats ()
    {
	local output="$1"
	shift
	local inputs=("$@")
	vcf_input_stream "${inputs[@]}" | "${VCF2STATS_BIN}" > "$output"
    }

    export -f vcf_input_stream prefix_log vcf_stats
    export VCF2STATS_BIN
    for ((i=0; i<${#INPUTS[@]}; i++)); do
	printf "set -exo pipefail; vcf_stats %q %q |& prefix_log vcfstats.%05d\n" "${STATS_FILES[i]}" "${INPUTS[i]}" "$i"
    done | "${PARALLEL_BIN}" --line-buffer

    echo "Columns extracted. Excerpt 0:" >&2
    head -n 10 "${STATS_FILES[0]}" >&2
    echo "Excerpt last (${#STATS_FILES[@]}):" >&2
    head -n 10 "${STATS_FILES[-1]}" >&2
)

(
    :
    : Calculate the distribution of DP values
    :
    set -eo pipefail

    {
	{
	    # include header from first file
	    cat "${STATS_FILES[0]}"
	    for ((i=1; i<${#STATS_FILES[@]}; i++)); do
		# skip header
		tail -n +2 "${STATS_FILES[i]}"
	    done
	} | {
	    calculate_distribution "DP" > "${workdir}/distribution.DP.txt";
	}
    } |& {
	prefix_log calculate_distribution yes
    }
)

echo "DP distribution:"
cat "${workdir}/distribution.DP.txt"
read N SUM_DP SUM_SQ_DP MEAN_DP POP_VAR POP_STDDEV SAM_VARIANCE SD_DP REST < <(tail -n1 "${workdir}/distribution.DP.txt")

(
    :
    : Apply filter to stats files.
    :
    # The output is CHR POS for each site that passes the filter.
    #
    function apply_filter ()
    {
	local stats_file="$1"
	local out_file="$2"
	shift 2
	local filter_args=("$@")
	wc -l "${stats_file}"
	"${AWK_BIN}" -f "${workdir}/filter.awk" "${filter_args[@]}" "${stats_file}" > "${out_file}" || {
	    echo "Failed to filter ${stats_file}" >&2
	    return 2
	}

	:
	: Sample output
	:
	wc -l "${out_file}"
	head "${out_file}"
    }
    export workdir AWK_BIN
    export -f apply_filter prefix_log
    generate_filter_prog "${PARAM_FILE}" > "${workdir}/filter.awk"
    :
    : Generated filter:
    :
    "${AWK_BIN}" '{printf("%03d %s\n", NR, $0)}' < "${workdir}/filter.awk"

    for ((i=0; i<${#STATS_FILES[@]}; i++)); do
	printf "set -exo pipefail; apply_filter %q %q -vMEAN_DP=%s -vSD_DP=%s |& prefix_log apply_filter.%05d yes\n" \
	       "${STATS_FILES[i]}" "${PASS_FILES[i]}" "${MEAN_DP}" "${SD_DP}" "$i"
    done | "${PARALLEL_BIN}" --line-buffer
)

(
    :
    : Create VCF files with selected sites.
    :

    # The PASS_FILES contain a list of "chr pos type" tuples for sites that
    # have passed the filter. We make new files containing just those
    # sites.
    function select_sites ()
    {
	local snpout="$1" indelout="$2"
	shift 2
	local sites=("$@")
	"${SELECTSITES_BIN}" >("${BGZIP_BIN}" -c > "$snpout") >("${BGZIP_BIN}" -c > "$indelout") "${sites[@]}"
    }
    export -f vcf_input_stream select_sites prefix_log
    export SELECTSITES_BIN BGZIP_BIN
    for ((i=0; i<${#INPUTS[@]}; i++)); do
	printf "set -exo pipefail; ( vcf_input_stream %q | select_sites %q %q %q; ) |& prefix_log selectsites.%05d yes\n" \
	       "${INPUTS[i]}" "${GOLD_SNPS[i]}" "${GOLD_INDELS[i]}" "${PASS_FILES[i]}" "$i"
    done | "${PARALLEL_BIN}" --line-buffer
)

(
    {
	# build snps gather argument
	if [[ -n "$REFERENCE" ]]; then
	    echo "-R $REFERENCE"
	fi
	echo "--TMP_DIR ${workdir}/gatk"
	for ((i=0; i<${#GOLD_SNPS[@]}; i++)); do
	    echo "-I ${GOLD_SNPS[i]}"
	done
    } > "${workdir}/snps.gather.arguments"

    {
	# build indels gather argument
	if [[ -n "$REFERENCE" ]]; then
	    echo "-R $REFERENCE"
	fi
	echo "--TMP_DIR ${workdir}/gatk"
	for ((i=0; i<${#GOLD_INDELS[@]}; i++)); do
	    echo "-I ${GOLD_INDELS[i]}"
	done
    } > "${workdir}/indels.gather.arguments"

    java_opts="--java-options \"-Xmx$((MAX_MEM_MB/2))m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true\""

    export -f prefix_log
    {
	echo "set -exo pipefail; HOME=${workdir}/gatk /gatk/gatk GatherVcfs ${java_opts} --arguments_file ${workdir}/snps.gather.arguments" \
	     "-O ${workdir}/gold.snps.vcf.gz |& prefix_log snps.gather"
	echo "set -exo pipefail; HOME=${workdir}/gatk /gatk/gatk GatherVcfs ${java_opts} --arguments_file ${workdir}/indels.gather.arguments" \
	     "-O ${workdir}/gold.indels.vcf.gz |& prefix_log indels.gather"
    } | "${PARALLEL_BIN}" --line-buffer

    {
	echo "set -exo pipefail; HOME=${workdir}/gatk /gatk/gatk IndexFeatureFile ${java_opts} -F ${workdir}/gold.snps.vcf.gz |& prefix_log snps.index"
	echo "set -exo pipefail; HOME=${workdir}/gatk /gatk/gatk IndexFeatureFile ${java_opts} -F ${workdir}/gold.indels.vcf.gz |& prefix_log indels.index"
    } | "${PARALLEL_BIN}" --line-buffer

    if [[ -n "$REFERENCE" ]]; then
	{
	    echo "set -exo pipefail; HOME=${workdir}/gatk /gatk/gatk ValidateVariants ${java_opts} -R ${REFERENCE} --validation-type-to-exclude IDS -V ${workdir}/gold.snps.vcf.gz |& prefix_log snps.validate"
	    echo "set -exo pipefail; HOME=${workdir}/gatk /gatk/gatk ValidateVariants ${java_opts} -R ${REFERENCE} --validation-type-to-exclude IDS -V ${workdir}/gold.indels.vcf.gz |& prefix_log indels.validate"
	} | "${PARALLEL_BIN}" --line-buffer
    fi
)

find "${workdir}"
SNPTMP=$(mktemp -p "$(dirname "${OUTPUTS[0]}")" tmp.midas.snp.XXXXXX)
INDELTMP=$(mktemp -p "$(dirname "${OUTPUTS[1]}")" tmp.midas.indel.XXXXXX)
TMPFILES=( "$SNPTMP"{,.tbi} "$INDELTMP"{,.tbi} )

# cp, possibly across filesystems
mv -v -- "${workdir}"/gold.snps.vcf.gz  "${SNPTMP}"
mv -v -- "${workdir}"/gold.indels.vcf.gz "${INDELTMP}"
mv -v -- "${workdir}"/gold.snps.vcf.gz.tbi "${SNPTMP}".tbi
mv -v -- "${workdir}"/gold.indels.vcf.gz.tbi "${INDELTMP}".tbi

# local. quick
mv -v -- "${SNPTMP}"       "${OUTPUTS[0]}"
mv -v -- "${SNPTMP}".tbi   "${OUTPUTS[0]}".tbi
mv -v -- "${INDELTMP}"     "${OUTPUTS[1]}"
mv -v -- "${INDELTMP}".tbi "${OUTPUTS[1]}".tbi
TMPFILES=()

END_TS=$(printf "%(%s)T" -1)

echo "Produced gold set in $(((END_TS - START_TS + 59)/60)) minutes"
