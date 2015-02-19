#!/bin/sh
# umask 000
for i in `seq 8`; do
    echo starting worker $i;
    ./crispor.cgi --worker >> temp/worker$i.log 2>&1 &
done
