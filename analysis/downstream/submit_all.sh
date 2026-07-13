#!/bin/bash
cd "$(dirname "$0")"
mkdir -p results
for job in bins quality leiden proportions; do
    bash submit.sh "$job"
done
