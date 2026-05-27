#!/bin/bash
cd "$(dirname "$0")"
shopt -s nullglob
for pidfile in log/sub_*.pid; do
    pid=$(cat "$pidfile")
    name=$(basename "$pidfile" .pid)
    if kill "$pid" 2>/dev/null; then
        echo "stopped $name (PID $pid)"
    else
        echo "$name (PID $pid) not running"
    fi
    rm -f "$pidfile"
done
