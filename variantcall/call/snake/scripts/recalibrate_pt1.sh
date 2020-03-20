#!/bin/bash

#
# Implements model generation for variant recalibration.
# This does a single read pass over all raw VCFs
#
# https://software.broadinstitute.org/gatk/documentation/article.php?id=2805
#
# This needs to be run from a GATK singularity container.
#
# It produces a .tgz file usable as the input to recalibrate_pt2.sh
#
# usage:
#
#  recalibrate_pt1.sh [--tmpdir TMPDIR] --gold-snps GOLDSNPS.VCF.GZ --gold-indels GOLDINDELS.VCF.GZ OUTPUT INPUTVCF+
#
#  If INPUTVCF starts with "@", then the parameter is assumed to be a
#  filename containing a list of vcf input filenames.
set -eo pipefail

CFGTMP=/tmp
INPUTS=()
OUTPUT=""
MAX_MEM_MB=$((16*1024))
GOLD_SNPS=""
GOLD_INDELS=""

# shut up perl.
#
# the LANG environment variable is propagated from the login node into
# the compute job. but a different set of locales are installed.
export LANG=C

: environment checks
hostname -f

#
# GATK needs this for plotting recalibration tables
#
if [[ -n "${CONDA_PREFIX:-}" ]]; then
    if [[ -n "${PATH:-}" ]]; then
	export PATH="${CONDA_PREFIX}/bin:$PATH"
    else
	export PATH="${CONDA_PREFIX}/bin"
    fi
fi

# Force early fail if the library is missing.  The GATK container
# includes R, but not necessary libraries, unfortunately. We provide
# them via conda.
which Rscript
which R
Rscript -e "library(ggplot2)"

echo process $$

while [[ $# -gt 0 ]]; do
    arg="$1"
    shift
    case "$arg" in
	--gold-snps)
	    GOLD_SNPS="$1"
	    shift
	    ;;
	--gold-indels)
	    GOLD_INDELS="$1"
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

set -x

if [[ -z "$OUTPUT" ]]; then
    echo "Missing output file" >&2
    exit 1
fi

if [[ "${#INPUTS[@]}" -lt 1 ]]; then
    echo "Missing input file(s)" >&2
    exit 1
fi

if [[ -n "${GOLD_SNPS}" && ! -e "${GOLD_SNPS}" ]]; then
    echo "Invalid SNPS gold set." >&2
    exit 1
fi

if [[ -n "${GOLD_INDELS}" && ! -e "${GOLD_INDELS}" ]]; then
    echo "Invalid INDELS gold set." >&2
    exit 1
fi

if [[ -z "${GOLD_SNPS}" && -z "${GOLD_INDELS}" ]]; then
    echo "Must provide at least one gold set" >&2
    exit 1
fi

if [[ -n "${GOLD_INDELS}" && -n "${GOLD_SNPS}" ]]; then
    # Running two of the same tool in parallel. And keep some extra.
    MAX_MEM_MB=$((MAX_MEM_MB / 2 - 2048))
else
    MAX_MEM_MB=$((MAX_MEM_MB - 1024))
fi

if [[ ${MAX_MEM_MB} -lt 1024 ]]; then
    echo "Provide more ram to --max-mem-mb <AMOUNT_MB>" >&2
    exit 1
fi

JAVA_OPTIONS="-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx${MAX_MEM_MB}m"


XTMP="${SLURM_TMPDIR:-$CFGTMP}" # use local disk if on cluster
mkdir -p "$XTMP" "$XTMP"/gatk
workdir=$(mktemp -d -p $XTMP tmp.recalibrate.XXXXXX);
trap 'rm --one-file-system -r ${workdir}' EXIT;

TMPSNPS="${workdir}/gold.snps.vcf.gz"
TMPINDELS="${workdir}/gold.indels.vcf.gz"
if [[ -n "${GOLD_SNPS}"   ]]; then cp "${GOLD_SNPS}"   "$TMPSNPS"; fi
if [[ -n "${GOLD_INDELS}" ]]; then cp "${GOLD_INDELS}" "$TMPINDELS"; fi

# Compute index files (.vcf.gz.tbi)
(
    [[ -n "${GOLD_SNPS}"   ]] && echo HOME="${workdir}" /gatk/gatk IndexFeatureFile -F "${TMPSNPS}"
    [[ -n "${GOLD_INDELS}" ]] && echo HOME="${workdir}" /gatk/gatk IndexFeatureFile -F "${TMPINDELS}" || :
) | "${CONDA_PREFIX}/bin/parallel" -j2 --line-buffer


argsfile="${workdir}/arguments.txt"
echo "building list of arguments" >&2
(
    set +x
    count=0
    for i in "${INPUTS[@]}"; do
	echo --variant "$i"
	let count=count+1
    done > "$argsfile"
    echo "$count vcfs added" >&2
)


resname=GOLD
resparams=known=false,training=true,truth=true,prior=10.0

# "/localscratch/jslegare.9353812.0/tmp.recalibrate.SLptUP/snp.recal.Rscript.pdf"
snp_rscript="snp.recal.Rscript"
snp_rscript_pdf="snp.recal.Rscript.pdf"
snp_recal="snp.recal.vcf.gz"
snp_recal_tbi="snp.recal.vcf.gz.tbi"
snp_model="snp.model.report"
snp_tranches="snp.tranches"
snp_tranches_pdf="snp.tranches.pdf"
indel_rscript="indel.recal.Rscript"
indel_rscript_pdf="indel.recal.Rscript.pdf"
indel_recal="indel.recal.vcf.gz"
indel_recal_tbi="indel.recal.vcf.gz.tbi"
indel_tranches="indel.tranches"
indel_tranches_pdf="indel.tranches.pdf"

# output
# snp.recal.Rscript
# snp.recal.Rscript.pdf
# snp.recal.vcf.gz
# snp.recal.vcf.gz.tbi
# snp.tranches
# snp.tranches.pdf

outputs=()
[[ -n "${GOLD_SNPS}"   ]] && outputs+=( "${snp_rscript}"   "${snp_rscript_pdf}"   "${snp_recal}"   "${snp_recal_tbi}"   "${snp_tranches}" "${snp_tranches_pdf}"  )  #"${snp_model}"
[[ -n "${GOLD_INDELS}" ]] && outputs+=( "${indel_rscript}" "${indel_rscript_pdf}" "${indel_recal}" "${indel_recal_tbi}" "${indel_tranches}" )

function variant_recal ()
{
    local logprefix="$1"
    local LINE
    shift
    set -o pipefail
    HOME="${workdir}" /gatk/gatk VariantRecalibrator \
	--java-options "$JAVA_OPTIONS" \
	--TMP_DIR "${XTMP}/gatk" \
	"$@" |& ( set +x; while read LINE; do printf "%s%s\n" "$logprefix" "$LINE"; done )
}

#https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/5585cdf7877104f2c61b2720ddfe7235f2fad577/JointGenotypingWf.wdl#L489
# /usr/gitc/gatk --java-options "-Xmx100g -Xms100g" \
#    VariantRecalibrator \
#    -V ${sites_only_variant_filtered_vcf} \
#    -O ${recalibration_filename} \
#    --tranches-file ${tranches_filename} \
#    --trust-all-polymorphic \
#    -tranche ${sep=' -tranche ' recalibration_tranche_values} \
#    -an ${sep=' -an ' recalibration_annotation_values} \
#    -mode SNP \
#    --sample-every-Nth-variant ${downsampleFactor} \
#    --output-model ${model_report_filename} \
#    --max-gaussians 6 \
#    -resource hapmap,known=false,training=true,truth=true,prior=15:${hapmap_resource_vcf} \
#    -resource omni,known=false,training=true,truth=true,prior=12:${omni_resource_vcf} \
#    -resource 1000G,known=false,training=true,truth=false,prior=10:${one_thousand_genomes_resource_vcf} \
#    -resource dbsnp,known=true,training=false,truth=false,prior=7:${dbsnp_resource_vcf}

# INDELS
# /usr/gitc/gatk --java-options "-Xmx24g -Xms24g" \
# 	       VariantRecalibrator \
#       -V ${sites_only_variant_filtered_vcf} \
#       -O ${recalibration_filename} \
#       --tranches-file ${tranches_filename} \
#       --trust-all-polymorphic \
#       -tranche ${sep=' -tranche ' recalibration_tranche_values} \
#       -an ${sep=' -an ' recalibration_annotation_values} \
#       -mode INDEL \
#       --max-gaussians 4 \
#       -resource mills,known=false,training=true,truth=true,prior=12:${mills_resource_vcf} \
#       -resource axiomPoly,known=false,training=true,truth=false,prior=10:${axiomPoly_resource_vcf} \
#       -resource dbsnp,known=true,training=false,truth=false,prior=2:${dbsnp_resource_vcf}

export JAVA_OPTIONS
export workdir
export XTMP
export -f variant_recal

# RSCript invocations:
# Executing: Rscript /localscratch/jslegare.9353812.0/tmp.recalibrate.SLptUP/snp.recal.Rscript
# Executing: Rscript (resource)org/broadinstitute/hellbender/tools/walkers/vqsr/plot_Tranches.R /localscratch/jslegare.9353812.0/tmp.recalibrate.SLptUP/snp.tranches 2.15
(
    # STEP 1 - build model for snps
    [[ -n "${GOLD_SNPS}" ]] && \
	echo set -x";" variant_recal SNP_ --arguments_file "$argsfile" \
	     --resource "${resname},${resparams}:${TMPSNPS}" \
	     --mode SNP -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
	     --trust-all-polymorphic \
	     --truth-sensitivity-tranche 100.0 \
	     --truth-sensitivity-tranche  99.0 \
	     --truth-sensitivity-tranche  90.0 \
	     --truth-sensitivity-tranche  70.0 \
	     --truth-sensitivity-tranche  50.0 \
	     `: --output-model "$workdir/${snp_model}"` \
	     --max-gaussians 6 \
	     --rscript-file  "$workdir/${snp_rscript}" \
	     --tranches-file "$workdir/${snp_tranches}" \
	     --output        "$workdir/${snp_recal}" || :

    # STEP 2 - build model for indels
    [[ -n "${GOLD_INDELS}" ]] && \
	echo set -x";" variant_recal IND_ --arguments_file "$argsfile" \
	     --resource "${resname},${resparams}:${TMPINDELS}" \
	     --mode INDEL -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR `: -an DP` \
	     --trust-all-polymorphic \
	     --truth-sensitivity-tranche 100.0 \
	     --truth-sensitivity-tranche  99.0 \
	     --truth-sensitivity-tranche  90.0 \
	     --truth-sensitivity-tranche  70.0 \
	     --truth-sensitivity-tranche  50.0 \
	     --max-gaussians 4 \
	     --bad-lod-score-cutoff -5.0 `: default is -5` \
	     --rscript-file  "$workdir/${indel_rscript}" \
	     --tranches-file "$workdir/${indel_tranches}" \
	     --output        "$workdir/${indel_recal}" || :

) | "${CONDA_PREFIX}/bin/parallel" -j2 --line-buffer

ls -lh "${workdir}"

( cd "${workdir}"
  if [[ -f /gatk/src/main/resources/org/broadinstitute/hellbender/tools/walkers/vqsr/plot_Tranches.R ]]; then
      cp /gatk/src/main/resources/org/broadinstitute/hellbender/tools/walkers/vqsr/plot_Tranches.R ./
      outputs+=("plot_Tranches.R")
  fi
  tar -czvf model.tgz "${outputs[@]}"
)
mv "${workdir}/model.tgz" "$OUTPUT"
