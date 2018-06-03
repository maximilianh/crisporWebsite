#!/bin/bash
set -e 
workCount=${1:-2}
echo `pwd`
if [[ `pwd` == *"crisporBeta"* ]] ; then
    workCount=1
fi
echo $workCount
mkdir -p log
for i in `seq $workCount`; do
    log=log/worker$i.log
    echo "---- NEW LOG ----" >> $log
    date >> $log
    pwd=`pwd`
    nohup ./crispor.py --worker $pwd $i >> $log 2>&1
    echo $! > log/worker$i.pid
    echo worker $i started, PID $!, logfile is $log;
done
