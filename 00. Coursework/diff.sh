#! /bin/bash

for t in 100 200 300 400 500
do
./Test/diff-output ./initial_output/output.dat$t output.dat$t
done
