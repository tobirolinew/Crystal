#!/bin/sh
#SBATCH --chdir=../clus/NN2-SAPT
#SBATCH --time=10:00
#SBATCH --cpus-per-task=60
#SBATCH --mem=8GB
#SBATCH --output=slurm.log
#SBATCH --error=slurm.log

BASE=~/Crystal

cat > inp.dat << EOF
Generic parameters:

  cut = 0 
  contrac = 0
  sys_name = rigidMC
  nthreads = ${SLURM_CPUS_PER_TASK}
  plot = 0

System-dependent parameters:

    Geometry:

      cluster  = none 
      monomer1 = ${BASE}/xyz/NN2/SAPT/NN.xyz
      monomer2 = ${BASE}/xyz/NN2/SAPT/NN.xyz

    Customize:

      custom_drive = normal
      ang_step[deg] = 50.0
      rad_step[A] = 0.5

    Cost:

      cost_label = Energy
      cost_unit  = cm-1
      cost_spec  = ${BASE}/par/N2d-SAPT.f2b
      cost_drive = f2b_homo

    Task:

      task_drive = msopt
      stac_point_file = stac.xyz
      max_targets = none
      max_ener = none
      max_order = none
      skip_flat = 1
EOF

module purge
module load oneapi/2024.2.1
${BASE}/crystal.x inp.dat > out.dat
rm -f "slurm.log" "done.log" "fort.2"
