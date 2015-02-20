#!/bin/sh
for i in `seq 4`; do
    log=temp/worker$i.log
    ./crispor.cgi --user www-data --worker >> $log 2>&1
    echo worker $i started, logfile is $log;
done
