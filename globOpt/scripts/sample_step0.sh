#!/bin/bash
noise_arg="--noise 0.003"
i=1; # arg0 is not used
for var in "$@"
do
    echo "$i $var"
    if [ "$var" == "--noise" ]; then
        tmp=$(($i+1))
        noise_arg="--noise ${!tmp}"
    fi
    ((i++))
done
echo "noisearg: ${noise_arg}"
cmd="Release/bin/gurobi_opt --sample-input ${noise_arg} --img $1"
echo "executing: $cmd"
${cmd}
