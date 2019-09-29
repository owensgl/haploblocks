#!/bin/bash

#
# Implements model generation for variant recalibration.
# This does a read pass over input raw VCFs and applies
# a model (produced with VariantRecalibrator).
#
#
# https://software.broadinstitute.org/gatk/documentation/article.php?id=2805
#
# This needs to be run from a GATK singularity container.
#
#  recalibrate_pt2.sh [--tmpdir TMPDIR] --model-snps MODEL.TGZ --model-indels MODEL.TGZ \
#                     [--tranches-snps TRANCHES] [--tranches-indels] TRANCHES OUTPUT INPUTVCF+
#
#
#  if --tranches-snp or --tranches-indel are provided, the file
#  specified will take precendence over the one specified via
#  --model-snp, or --model-indel, respectively.
#
#  MODEL.TGZ is assumed to contain a snp.recal and snp.tranches file
#  for snps. For indels, indel.recal, indel.tranches
#
#  If INPUTVCF starts with "@", then the parameter is assumed to be a
#  filename containing a list of vcf input filenames.
#
set -exo pipefail

CFGTMP=/tmp
INPUTS=()
OUTPUT=""
MAX_MEM_MB=$((4*1024))
MODEL_SNP=""
MODEL_INDEL=""
TRANCHES_SNP=""
TRANCHES_INDEL=""

while [[ $# -gt 0 ]]; do
    arg="$1"
    shift
    case "$arg" in
	--model-snps)
	    MODEL_SNP="$1"
	    shift
	    ;;
	--model-indels)
	    MODEL_INDEL="$1"
	    shift
	    ;;
	--tranches-snps)
	    TRANCHES_SNP="$1"
	    shift
	    ;;
	--tranches-indels)
	    TRANCHES_INDEL="$1"
	    shift
	    ;;
	--tmpdir)
	    CFGTMP="$1"
	    shift
	    ;;
	--max-mem-mb)
	    MAX_MEM_MB="$1"
	    shift
	    ;;
	@*)
	    # read filenames from file
	    infile="${arg##@}"
	    while read inputvcf rest; do
		INPUTS+=( "$inputvcf" )
	    done < "$infile"
	    ;;
	*)
	    if [[ -z "$OUTPUT" ]]; then
		OUTPUT="$arg"
	    else
		INPUTS+=( "$arg" )
	    fi
	    ;;
    esac
done

JAVA_OPTIONS="-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx${MAX_MEM_MB}m"

if [[ -z "$OUTPUT" ]]; then
    echo "Missing output file" >&2
    exit 1
fi

if [[ "${#INPUTS[@]}" -lt 1 ]]; then
    echo "Missing input file(s)" >&2
    exit 1
fi

if [[ -z "$MODEL_SNP" ]]; then
    echo "Missing SNP variant model" >&2
    exit 1
fi

if [[ -z "$MODEL_INDEL" ]]; then
    echo "WARN: No indel model provided. indels will not be annotated." >&2
fi

# Use SLURM_TMPDIR, if on cluster, to store model.
XTMP="${SLURM_TMPDIR:-$CFGTMP}"
# Do not use SLURM_TMPDIR for the final output. We run out of space
# when jobs get co-scheduled.
OTMP="${CFGTMP}"

mkdir -p "$XTMP" "$XTMP"/gatk "$OTMP"

workdir=$(mktemp -d -p "$XTMP" tmp.recalibrate_pt2_model.XXXXXX);
workdir_output=$(mktemp -d -p "$OTMP" tmp.recalibrate_pt2_output.XXXXXX);

trap 'find "${workdir}" || : ; find "${workdir_output}" || : ; rm --one-file-system -r -- "${workdir}" "${workdir_output}"' EXIT;


function unpack_model () # tgz outputdir basename.vcf.gz
{
    local outdir="$2" tgz="$1" vcf="$3"
    mkdir -p "$outdir"
    tar -C "$outdir" -xvf "$tgz"
    if [[ ! -f "$vcf" ]]; then
	echo "File $vcf is missing." >&2
	return 1
    fi
    if [[ ! -f "$vcf".tbi ]]; then
	HOME="${workdir}" /gatk/gatk IndexFeatureFile -F "${vcf}"
    fi
}

# Example foo.tranches file content
#
# # Variant quality score tranches file
# # Version number 5
# targetTruthSensitivity,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,model,accessibleTruthSites,callsAtTruthSites,truthSensitivity
# 50.00,0,3299230,0.0000,2.2820,38.5404,VQSRTrancheSNP0.00to50.00,SNP,648748,324374,0.5000
# 70.00,0,18837779,0.0000,2.2139,20.0489,VQSRTrancheSNP50.00to70.00,SNP,648748,454123,0.7000
# 90.00,0,22418253,0.0000,2.1860,17.1817,VQSRTrancheSNP70.00to90.00,SNP,648748,583873,0.9000
# 99.00,0,47436477,0.0000,2.0826,10.9903,VQSRTrancheSNP90.00to99.00,SNP,648748,642260,0.9900
# 100.00,0,171274216,0.0000,1.7160,-39999.3077,VQSRTrancheSNP99.00to100.00,SNP,648748,648748,1.0000
function get_tranche_values () {
    # get a sorted list of tranche values recorded in the tranches file, in ascending order.
    local VAL REST
    local NR=0
    while IFS="," read VAL REST; do
	case "$VAL" in
	    "#"*)
		continue
		;;
	esac
	NR=$((NR + 1));
	if [[ "$NR" -eq 1 ]]; then
	    if [[ "$VAL" != "targetTruthSensitivity" ]]; then
		echo "tranches file format different than expected" >&2
		return 1
	    fi
	    continue
	fi
	echo "Found tranche $VAL" >&2
	echo "$VAL"
    done | python -c "import sys
floats = sorted([(float(line.strip()),line.strip()) for line in sys.stdin])
for _, strval in floats: print(strval)"
}


SNPDIR="${workdir}/model_snp"
INDELDIR="${workdir}/model_indel"

# Contained in the model.tgz
snp_rscript="${SNPDIR}/snp.recal.Rscript"
snp_recal="${SNPDIR}/snp.recal.vcf.gz"
snp_tranches="${SNPDIR}/snp.tranches"
snp_tranche_values=()
indel_rscript="${INDELDIR}/indel.recal.Rscript"
indel_recal="${INDELDIR}/indel.recal.vcf.gz"
indel_tranches="${INDELDIR}/indel.tranches"
indel_tranche_values=()

mkdir -p "${SNPDIR}" "${INDELDIR}"

if [[ ! -z "$MODEL_SNP" ]]; then
    unpack_model "${MODEL_SNP}" "${SNPDIR}" "${snp_recal}"
fi

if [[ ! -z "$MODEL_INDEL" ]]; then
    unpack_model "${MODEL_INDEL}" "${INDELDIR}" "${indel_recal}"
fi

argsfile="${workdir}/arguments.txt"
(
    for i in "${INPUTS[@]}"; do
	echo --variant "$i"
    done
) > "$argsfile"


if [[ -n "$TRANCHES_SNP" ]]; then
    # cli override
    snp_tranches="${TRANCHES_SNP}"
fi

if [[ -n "$TRANCHES_INDEL" ]]; then
    # cli override
    indel_tranches="${TRANCHES_INDEL}"
fi

if [[ -n "${MODEL_SNP}" ]]; then
    :
    : SNP tranches file
    :
    cat "${snp_tranches}"
    snp_tranche_values=( $(get_tranche_values < "${snp_tranches}") )
fi

if [[ -n "${MODEL_INDEL}" ]]; then
    :
    : INDEL tranches file
    :
    cat "${indel_tranches}"
    indel_tranche_values=( $(get_tranche_values < "${indel_tranches}") )
fi

previnput=( --arguments_file "$argsfile" )
nextoutput="${workdir_output}/tmp_tranches.pt1.vcf.gz"

:
: annotate SNPs with PASS or FILTER for the given sensitivity levels.
:
HOME="${workdir}" /gatk/gatk ApplyVQSR \
    --java-options "$JAVA_OPTIONS" \
    --TMP_DIR "${XTMP}/gatk" \
    "${previnput[@]}" \
    --mode SNP --truth-sensitivity-filter-level "${snp_tranche_values[0]}" \
    --tranches-file "${snp_tranches}"  \
    --recal-file "${snp_recal}" \
    --create-output-variant-index true \
    --output "${nextoutput}"

if [[ ! -z "$MODEL_INDEL" ]]; then
    previnput=( --variant "${nextoutput}" )
    nextoutput="${workdir_output}/tmp_tranches.pt2.vcf.gz"

    :
    : annotate INDELs with PASS or FILTER for the given sensitivity levels.
    :
    HOME="${workdir}" /gatk/gatk ApplyVQSR \
	--java-options "$JAVA_OPTIONS" \
	--TMP_DIR "${XTMP}/gatk" \
	"${previnput[@]}" \
	--mode INDEL --truth-sensitivity-filter-level "${indel_tranche_values[0]}" \
	--tranches-file "${indel_tranches}" \
	--recal-file "${indel_recal}" \
	--create-output-variant-index true \
	--output "${nextoutput}"
fi

mv "${nextoutput}" "$OUTPUT"
if [[ -e "${nextoutput}.tbi" ]]; then
    mv "${nextoutput}.tbi" "$OUTPUT".tbi
fi
