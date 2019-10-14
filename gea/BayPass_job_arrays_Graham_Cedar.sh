#!/bin/bash
# ---------------------------------------------------------------------
# SLURM script for BayPass (IS_STD model) as job arrays
# ---------------------------------------------------------------------
#SBATCH --account=def-yeaman
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=25000mb
#SBATCH --tasks=1
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shaghayegh.soudi@gmail.com

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
# ---------------------------------------------------------------------
echo ""
echo "Job Array ID / Job ID: $SLURM_ARRAY_JOB_ID / $SLURM_JOB_ID"
echo "This is job $SLURM_ARRAY_TASK_ID out of $SLURM_ARRAY_TASK_COUNT jobs."
echo ""

## 28.12.2018, script to run BayPass covariate model on Cedar 
## load necessary modules to run BayPass

module load nixpkgs/16.09  
module load intel/2016.4 
module load baypass/2.1


## Directories:
# Locate the input data
ROOT_DIR=/project/6004169/shsoudi/H_petiolaris_HA412/pet.can
GENO_DIR=${ROOT_DIR}/5K_genotype_files_pet_can
ENV_DIR=${ROOT_DIR}/pet.can.env.file
MATRIX_DIR=${ROOT_DIR}/pet.can.matrix.file


# Output directory
OUTPUT_DIR=${ROOT_DIR}/pet.can.Baypass.IS.STD.output.file

## Input file
BASENAME=baypass_geno_pet.can.HA412
GENO=${GENO_DIR}/${BASENAME}_${SLURM_ARRAY_TASK_ID}
ENV=${ENV_DIR}/data_climate_NA_H_petiolaris_CAN_2018_formatted_baypass.txt
MAT=${MATRIX_DIR}/pet.can.random1.anacore_mat_omega.out
 
 
i_baypass -npop 9 -gfile ${GENO} -efile ${ENV} -omegafile ${MAT} -nthreads $SLURM_CPUS_PER_TASK -outprefix ${OUTPUT_DIR}/petiolaris_CAN_anacovis2_${SLURM_ARRAY_TASK_ID}

 
 
# ---------------------------------------------------------------------
echo "Job finished with exit code $? at: `date`"
# ---------------------------------------------------------------------



