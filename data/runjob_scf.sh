#!/usr/bin/env sh
#
# runjob_scf.sh
#

alpha=(1e00 1e-1 1e-2)
beta=(1e00 1e-1 1e-2)
z=$2
source=$PATH_TO_1D2E_TESTS
echo "Atomic Number: "$z >> job.log
echo "Running SCF jobs in "$pwd >> job.log
echo "No. of threads being used: "$OMP_NUM_THREADS >> job.log

# DVR grid parameters
x0="15.0"
ngrid=(32 64 128 256 512 1024)
# ngrid2=(2048 4096 8192)

for a in ${alpha[@]};do
    for b in ${beta[@]};do
        cd "alpha_"$a"_beta_"$b;
        num_alpha=$(printf "%.8f\n" "$a")
        num_beta=$(printf "%.8f\n" "$b")
        num_z=$(printf "%.8f\n" "$z")
        for n in ${ngrid[@]};do
            cd "n_"$n;
            echo "Running SCF calculations in alpha_"$a"_beta_"$b"/n_"$n >> ../../job.log
            echo "Running Hartree-Fock" >> run.log
            echo $source/test_hartreefock -n $n -x0 $x0 -z $num_z -alpha $num_alpha -beta $num_beta >> run.log
            $source/test_hartreefock -n $n -x0 $x0 -z $num_z -alpha $num_alpha -beta $num_beta &> scf.log
            cd ..
        done
        cd ..
    done
done

