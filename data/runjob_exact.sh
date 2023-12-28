#!/usr/bin/env sh
#
# runjob_exact.sh
#

alpha=(1e00 1e-1 1e-2)
beta=(1e00 1e-1 1e-2)
z=$2
source=$PATH_TO_1D2E_TESTS
echo "Atomic Number: "$z >> job.log
echo "Running Exact Diagonalization jobs in "$pwd >> job.log
echo "No. of threads being used: "$OMP_NUM_THREADS >> job.log

# DVR grid parameters
x0="15.0"
ngrid=(32 64 128)

for a in ${alpha[@]};do
    for b in ${beta[@]};do
        cd "alpha_"$a"_beta_"$b;
        num_alpha=$(printf "%.8f\n" "$a")
        num_beta=$(printf "%.8f\n" "$b")
        num_z=$(printf "%.8f\n" "$z")
        for n in ${ngrid[@]};do
            cd "n_"$n;
            echo "Running Exact Diagonalization" >> run.log
            echo $source/test_2e_dvr -n $n -x0 $x0 -z $num_z -alpha $num_alpha -beta $num_beta >> run.log
            $source/test_2e_dvr -n $n -x0 $x0 -z $num_z -alpha $num_alpha -beta $num_beta &> exact.log

            echo "Computing the Wigner Intracule" >> run.log
            $source/test_intracules -i psi_exact.wfn -x0 $x0 -o wigner_intracule.wfn &> intracules.log
            echo $source/test_intracules -i psi_exact.wfn -x0 $x0 -o wigner_intracule.wfn &> intracules.log
            cd ..
        done
        cd ..
    done
done
