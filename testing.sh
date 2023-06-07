#!/bin/bash
echo "Start testing"

for i in {2..7}
do
    echo "Test $i"
    # Calculate 2 ^ i using bash
    POWER=$((2 ** $i))
    sage run.py -g $POWER
    sage run.py -t msidh -r 10 -f ./MSIDH_AES-$POWER.pickle
done