#
#
# Example with hello world 
# We are planning to add it as EvaluationSLURM() using fork() - exec()?
#
#

#!/bin/sh
#SBATCH --partition=general-compute
#SBATCH --time=00:30:00
#SBATCH --nodes=40
#SBATCH --ntasks-per-node=8
#SBATCH --constraint=IB
#SBATCH --mem=23000
# Memory per node specification is in MB. It is optional. 
# The default limit is 3000MB per core.
#SBATCH --job-name="hello_test"
#SBATCH --output=test-srun.out
#SBATCH --mail-user=oksana.shadura@cern.ch
#SBATCH --mail-type=ALL
##SBATCH --requeue
#Specifies that the job will be requeued after a node failure.
#The default is that the job will not be requeued.

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "Working directory = "$SLURM_SUBMIT_DIR

# Check  which modules to upload
module load intel/13.1
module load intel-mpi/4.1.3
module list
ulimit -s unlimited
#
# The initial srun will trigger the SLURM prologue on the compute nodes.
NPROCS=`srun --nodes=${SLURM_NNODES} bash -c 'hostname' |wc -l`
echo NPROCS=$NPROCS
echo "Launch helloworld with srun"
#The PMI library is necessary for srun
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
srun ./helloworld
#
echo "All Done!"