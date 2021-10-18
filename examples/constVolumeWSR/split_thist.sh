#!/bin/bash

thist_file="$1"

if [ ${thist_file}"X" == "X" ]
then
  echo "Usage: ${0} <thist_file>"
  exit 1
fi

nruns=`grep "...." ${thist_file} | tail -n 1 | awk '{print $1;}'`
for i in `seq 0 $nruns`
do
  grep "^\ *${i}\ \ *" ${thist_file} > run_id${i}.thist
done

