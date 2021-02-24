#!/bin/bash

# ezq.sh -- easily add pw.x runs to the job queue
# John Wilkinson 24/2/21. Based on the input file by Yinan Chen.

function print_help {
	cat < /dev/tty << EOF
EzQ (EasyQueue) -- Queue runs for pw.x on Redwood nodes with ease!
-----------------------------------------------------------------
Usage: ezq -m [memory in GB] -n [number of nodes] -p [number of thread pools] -t [number of threads] PWI_FILE

-m --memory: Memory required for the run. Default is 3 GB, max for Redwood is 6. You might get a JobHeldUser error in squeue if this is too high.

-n --nodes: Number of nodes requested. There are currently 19 on Redwood.

-p --pools: Number of thread pools (the argument -npool in pw.x). Default is 2; you can probably leave this alone.

-t --threads: Number of threads. Sets the OMP_NUM_THREADS variable, default is 1. 

PWI_FILE: pw.x input file. The output of this file will be saved as this file name with the extension .pwo
EOF
}

# set the variables of the run 
nproc=32 #number of processors per node  !!!! This will change with the new nodes

nodes=2  #total number of processors divided by number of processors per node
nthreads=1 #number of threads
npool=2 # number of pools

queue=redwood # queue type: short, long, bigmem   !!!! This will change with the new nodes

memory=3 # required memory in GB

# load in the ezq configuration (which gives the paths we want)
. $(dirname `which $0`)/ezq.conf

# parse the arguments
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -h|--help)
	    print_help
	    exit
	    ;;
    -m|--memory)
    memory="$2"
    shift # past argument
    shift # past value
    ;;
    -n|--nodes)
    nodes="$2"
    shift # past argument
    shift # past value
    ;;
    -p|--pools)
    npool="$2"
    shift # past argument
    shift # past value
    ;;
    -t|--threads)
    nthreads="$2" 
    shift # past argument
    shift # past value
    ;;
    *)    #  file name
    pwi_file="$1" # pwi file 
    shift # past argument
    ;;
esac
done

if [ -z "$pwi_file" ]
then
	echo "No pwi file found. Aborting";
	exit
fi


# now adjust the number of processors to match the number of nodes*proc per node
ntot=$(( nproc*nodes ))  #total number of processors

file_root="${pwi_file%.*}"
pwo_file="${file_root}.pwo"
run_file="${file_root}.run"

cat > /dev/tty << EOF
Setting a job up on Redwood with $nodes nodes, and hence a total of $ntot processors.
The memory per node is $memory GB, nthreads=$nthreads, npools=$npool
The pw.x input file is $pwi_file
The pw.x output will be saved to $pwo_file
EOF

cat > $run_file << EOF
#!/bin/bash

module load intel-compilers
module load openmpi
module load fftw/3.3.8-intel

export OMP_NUM_THREADS=$nthreads

$mpi -n $ntot -m cyclic --mpi=pmi2 $qe_path/pw.x -npool $npool < $pwi_file > $pwo_file
EOF

chmod a+x $run_file
addqueue -s -c "$file_root" -n "$nodes"x"$nproc" -m $memory -q $queue ./$run_file

