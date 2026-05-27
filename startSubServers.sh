#!/bin/bash
set -e
cd "$(dirname "$0")"
mkdir -p log

while IFS='|' read -r name port venv; do
    log=log/sub_${name}.log
    pid=log/sub_${name}.pid
    echo "---- NEW LOG ----" >> "$log"
    date >> "$log"
    nohup ./$venv/bin/python3 startSubServer.py "$port" "$name" >> "$log" 2>&1 &
    echo $! > "$pid"
    echo "subserver $name (port $port, venv $venv) started, PID $!, log $log"
done < <(python3 -c "
from subserversConf import SUBSERVERS
for n, c in SUBSERVERS.items():
    print(f\"{n}|{c['port']}|{c['venv']}\")
")
