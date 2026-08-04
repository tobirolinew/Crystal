#!/bin/sh
#SBATCH --chdir=../clus/LJ5
#SBATCH --time=10:00
#SBATCH --cpus-per-task=60
#SBATCH --mem=10GB
#SBATCH --output=slurm.log
#SBATCH --error=slurm.log

BASE=~/Crystal

cat > inp.dat << EOF
Generic parameters:

  cut = 0.0 
  contrac = 0.01
  sys_name = rigidMC
  nthreads = ${SLURM_CPUS_PER_TASK}
  plot = 0

System-dependent parameters:

    Geometry:

      cluster  = ${BASE}/xyz/LJn/LJ5.xyz
      monomer1 = ${BASE}/xyz/LJn/LJ.xyz
      monomer2 = ${BASE}/xyz/LJn/LJ.xyz
      monomer3 = ${BASE}/xyz/LJn/LJ.xyz
      monomer4 = ${BASE}/xyz/LJn/LJ.xyz
      monomer5 = ${BASE}/xyz/LJn/LJ.xyz

    Customize:

      custom_drive = normal
      ang_step[deg] = 45.0
      rad_step[A] = relax

    Cost:

      cost_label = Energy
      cost_unit  = cm-1
      cost_spec  = ${BASE}/par/LJ2.f2b
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
