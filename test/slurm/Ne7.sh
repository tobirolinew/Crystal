#!/bin/sh
#SBATCH --chdir=../clus/Ne7
#SBATCH --time=10:00
#SBATCH --cpus-per-task=60
#SBATCH --mem=8GB
#SBATCH --output=slurm.log
#SBATCH --error=slurm.log

BASE=~/Crystal

cat > inp.dat << EOF
Generic parameters:

  cut = -369.6 
  contrac = 0.9
  sys_name = rigidMC
  nthreads = ${SLURM_CPUS_PER_TASK}
  plot = 0

System-dependent parameters:

    Geometry:

      cluster  = ${BASE}/xyz/Nen/Ne7_min2.xyz
      monomer1 = ${BASE}/xyz/Nen/Ne.xyz
      monomer2 = ${BASE}/xyz/Nen/Ne.xyz
      monomer3 = ${BASE}/xyz/Nen/Ne.xyz
      monomer4 = ${BASE}/xyz/Nen/Ne.xyz
      monomer5 = ${BASE}/xyz/Nen/Ne.xyz
      monomer6 = ${BASE}/xyz/Nen/Ne.xyz
      monomer7 = ${BASE}/xyz/Nen/Ne.xyz

    Customize:

      custom_drive = normal
      ang_step[deg] = 90.0
      rad_step[A] = relax

    Cost:

      cost_label = Energy
      cost_unit  = cm-1
      cost_spec  = ${BASE}/par/Ne2.f2b
      cost_drive = f2b_homo

    Task:

      task_drive = msopt
      stac_point_file = stac.xyz
      max_targets = none
      max_ener = none
      max_order = 0
      skip_flat = 0
EOF

module purge
module load oneapi/2024.2.1
${BASE}/crystal.x inp.dat > out.dat
rm -f "slurm.log" "done.log" "fort.2"
