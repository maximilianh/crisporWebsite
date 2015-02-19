#!/bin/sh
for i in `seq 4`; do
    echo starting worker $i;
    ./crispor.cgi --worker >> temp/worker$i.log 2>&1 &
done
