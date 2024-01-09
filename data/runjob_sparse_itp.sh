#!/usr/bin/env sh
#
# runjob_sparse_itp.sh
#

alpha=(1e00 1e-1 1e-2)
beta=(1e00 1e-1 1e-2)
z=$2
source=$PATH_TO_1D2E_TESTS
echo "Atomic Number: "$z >> job.log
echo "Running Sparse ITP jobs in "$pwd >> job.log
echo "No. of threads being used: "$OMP_NUM_THREADS >> job.log

# DVR grid parameters
x0="15.0"
ngrid=(256) #(512 1024)

# time steps for itp
short_tstep=0.0001
etol=1E-10

for a in ${alpha[@]};do
    for b in ${alpha[@]};do
        cd "alpha_"$a"_beta_"$b;
        num_alpha=$(printf "%.8f\n" "$a")
        num_beta=$(printf "%.8f\n" "$b")
        num_z=$(printf "%.8f\n" "$z")
        for n in ${ngrid[@]};do
            cd "n_"$n;
            # echo "Running Sparse-ITP calculations in alpha_"$a"_beta_"$b"/n_"$n >> ../../job.log
            # echo "Running Sparse-ITP" >> run.log
            # $source/test_sparse_itp -n $n -x0 $x0 -z $num_z -alpha $num_alpha -beta $num_beta -dt $short_tstep -etol $etol -print_nstep 1000 &> sparse_itp.log
            # echo $source/test_sparse_itp -n $n -x0 $x0 -z $num_z -alpha $num_alpha -beta $num_beta -dt $short_tstep -etol $etol -print_nstep 1000 &> sparse_itp.log
		                      
            echo "Computing the Wigner Intracule" >> run.log
            $source/test_intracules -i psi_itp.wfn -x0 $x0 -o itp_wigner_intracule.wfn &> intracules.log
            echo $source/test_intracules -i psi_itp.wfn -x0 $x0 -o itp_wigner_intracule.wfn &> intracules.log
            cd ..
        done
        cd ..
    done
done
