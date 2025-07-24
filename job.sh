#!/bin/sh

ifort -c bacs.f90
ifort -c functions.f90
ifort -c ranpack.f90
ifort -c Hamiltonian.f90

ifort -lpthread -lm -qmkl=sequential -Wl bacs.o functions.o ranpack.o Hamiltonian.o bdg_real.f90


d_delta = 0.01
for iSeed in $(seq 101 1 105)
for d in $(seq 0 1 5); do
    jobname="job_iSeed${iSeed}_D${ds}"

    cat > ${jobname}.sh <<EOF
#!/bin/sh`
#SBATCH --partition Inclusive
#SBATCH --nodelist Ubuntu
#SBATCH --cpus-per-task 1
#SBATCH --job-name=${jobname}
#SBATCH --time infinite
#SBATCH --output out_dir/MPI_%j.out

mkdir iSeed${iSeed}_d${ds}

./a.out $iSeed $ds $d_delta
 make clean

EOF

    sbatch ${jobname}.sh
done
done

