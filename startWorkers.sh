#!/bin/sh
set -e 
if [ "$EUID" -ne 0 ]
  then echo "must be run with sudo"
  exit
fi
workCount=4
if [ `pwd` == "/var/www/crisporTest" -o `pwd` == "/var/www/crisporMax" ] ; then
    workCount=1
fi
echo $workCount
for i in `seq $workCount`; do
    log=log/worker$i.log
    date > $log
    ./crispor.cgi --user www-data --worker `pwd` $i >> $log 2>&1
    echo worker $i started, PID $!, logfile is $log;
done
