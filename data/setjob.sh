#!/usr/bin/env sh
#
# setjob.sh
#

alpha=(1e00 1e-1 1e-2)
beta=(1e00 1e-1 1e-2)
ngrid=(32 64 128 256 512 1024 2048 4096 8192)

for a in ${alpha[@]};do
    for b in ${beta[@]};do
        mkdir "alpha_"$a"_beta_"$b;
        cd "alpha_"$a"_beta_"$b;
        for n in ${ngrid[@]};do
            mkdir "n_"$n;
            cd "n_"$n;
            touch run.log
            cd ..
        done
        cd ..
    done
done
