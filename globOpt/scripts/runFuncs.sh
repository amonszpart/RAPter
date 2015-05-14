#!/bin/bash

# save run.log arguments to "run.log"
function save_args() {
	logfile="run.log"
	args=("$@") 
	# get number of elements 
	ELEMENTS=${#args[@]} 

	echo -n "[$(date +%D\ %T)] " >> $logfile
	# echo each element in array  
	# for loop 
	for (( i=0;i<$ELEMENTS;i++)); do 
	    echo -n "${args[${i}]} " >> $logfile
	done
	# endline
	echo -e -n "\n" >> $logfile
}

function countLines()
{
    echo `wc -l $1 | cut -f1 -d' '`;
}

# only works for |ext|==3
function cutExt()
{
	echo `echo $1 | rev | cut -c 5- | rev`
}

# show command before run
# stop script when the command fails
function my_exec() {
	echo "__________________________________________________________";
	echo -e "\n\n[CALLING] $1";
    echo -e $1"\n" >>$execLog
    if ! $dryRun ; then
        eval $1;
    fi
  	if [ "$?" -ne "0" ]; then
	    echo "Error detected ($?). ABORT."
	    exit 1
  	fi
}

# won't halt, just report return value
function my_exec2() {
        echo "__________________________________________________________";
        echo -e "\n\n[CALLING] $1";
        echo -e $1"\n" >>$execLog
        if ! $dryRun ; then
            eval $1;
        fi
        ret=$?;
        if [ "$ret" -ne "0" ]; then
            echo "Call returned $ret"
        fi

        myresult=$ret;
}

function my_mult() {
	#eval "bc <<< 'scale=5; $1 * $2'"
	eval "bc <<< 'scale=5; $1 * $2'" | awk '{printf "%.6f", $0}'
}