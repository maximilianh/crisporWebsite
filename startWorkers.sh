#!/bin/sh
if [ "$EUID" -ne 0 ]
  then echo "must be run with sudo"
  exit
fi
for i in `seq 4`; do
    log=log/worker$i.log
    ./crispor.cgi --user www-data --worker `pwd` > $log 2>&1
    echo worker $i started, logfile is $log;
done
