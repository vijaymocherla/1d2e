#!/usr/bin/env sh
#
# runjob_atoms.sh
# -------------------------------------------------------
#  
#
# 

alpha=(1e00 1e-1 1e-2)
beta=(1e00 1e-1 1e-2)
atom_nos=(1 ) #2 3 4 5 6)
source="/home/sai/Desktop/1d2e/build/tests"
# DVR grid parameters
x0="15.0"
small_n=(32 64 128)
large_n=(256 512 1024)

# time steps for itp
long_tstep=0.001 
short_tstep=0.0001
etol=1E-12
print_nstep=100

for z in ${atom_nos[@]};do
    mkdir "z_"$z;
    cd "z_"$z;
    for a in ${alpha[@]};do
        for b in ${alpha[@]};do
            mkdir "alpha_"$a"_beta_"$b;
            cd "alpha_"$a"_beta_"$b;
            num_alpha=$(printf "%.8f\n" "$a")
            num_beta=$(printf "%.8f\n" "$b")
            num_z=$(printf "%.8f\n" "$z")
            for n in ${small_n[@]};do
                mkdir "n_"$n;
                cd "n_"$n;
                echo "Running Hartree-Fock" >> run.log
                $source/test_hartreefock -n $n -x0 $x0 -z $num_z -alpha $num_alpha -beta $num_beta &> scf.log
                echo $source/test_hartreefock -n $n -x0 $x0 -z $num_z -alpha $num_alpha -beta $num_beta &> scf.log

                echo "Running On-the-fly ITP" >> run.log                
                $source/test_2e_dvr -n $n -x0 $x0 -z $num_z -alpha $num_alpha -beta $num_beta -dt $long_tstep -etol $etol -print_nstep $print_nstep &> otf_itp.log
                echo $source/test_2e_dvr -n $n -x0 $x0 -z $num_z -alpha $num_alpha -beta $num_beta -dt $long_tstep -etol $etol -print_nstep $print_nstep &> otf_itp.log

                echo "Running Sparse-ITP" >> run.log
                $source/test_sparse_itp -n $n -x0 $x0 -z $num_z -alpha $num_alpha -beta $num_beta -dt $long_tstep -etol $etol -print_nstep $print_nstep &> sparse_itp.log
                echo $source/test_sparse_itp -n $n -x0 $x0 -z $num_z -alpha $num_alpha -beta $num_beta -dt $long_tstep -etol $etol -print_nstep $print_nstep &> sparse_itp.log

                echo "Computing the Wigner Intracule" >> run.log
                $source/test_intracules -i psi_itp.wfn -x0 $x0 -o wigner_intracule.wfn &> intracules.log
                echo $source/test_intracules -i psi_itp.wfn -x0 $x0 -o wigner_intracule.wfn &> intracules.log
                cd ..
            done
            cd ..
        done
    done
    cd ..
done