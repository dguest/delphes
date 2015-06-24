#!/usr/bin/env bash

set -eu

LOGFILE=memory-use.log
rm -f $LOGFILE

while true ; do
    # get time and memory use in KB
    ps aux | grep -v ' awk ' | awk '/dguest.*Delphes/ {print systime() " " $6}'
    sleep 1
done >> $LOGFILE

