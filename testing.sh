#!/bin/bash
echo "Start testing"

for i in {2..6}
do
    echo "Test $i"
    # Calculate 2 ^ i using bash
    POWER=$((2**$i))
    sage run.py -g $POWER
    # sage run.py -t msidh -r 3 -f ./MSIDH_AES-$i.pickle
done