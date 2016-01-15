#!/bin/bash
set -e 
workCount=${1:-4}
echo `pwd`
if [ `pwd` == "/var/www/crisporTest" -o `pwd` == "/var/www/crisporMax" -o `pwd` == "/home/www/crisporMax" ] ; then
    workCount=1
fi
echo $workCount
mkdir -p log
for i in `seq $workCount`; do
    log=log/worker$i.log
    date > $log
    ./crispor.cgi --worker `pwd` $i >> $log 2>&1
    echo worker $i started, PID $!, logfile is $log;
done
