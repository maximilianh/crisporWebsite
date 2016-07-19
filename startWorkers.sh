#!/bin/bash
set -e 
workCount=${1:-2}
echo `pwd`
if [[ `pwd` == "*crisporDev*" ]] ; then
    workCount=1
fi
echo $workCount
mkdir -p log
for i in `seq $workCount`; do
    log=log/worker$i.log
    date > $log
    ./crispor.py --worker `pwd` $i >> $log 2>&1
    echo worker $i started, PID $!, logfile is $log;
done
