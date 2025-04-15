#!/bin/bash
LOGFILE="./cpu_usage.log"

while true; do
    echo "==== $(date) ====" >> "${LOGFILE}"
    mpstat -P ALL 1 1 >> "${LOGFILE}"
    sleep 10
done
