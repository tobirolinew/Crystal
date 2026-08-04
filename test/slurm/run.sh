for f in *.sh
do
   if [[ $f == "run.sh" ]]; then continue; fi
   echo $f
   sbatch $f
done
