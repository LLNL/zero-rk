#!/bin/bash

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../
rm -rf outputs
mkdir outputs

setYMLScalar()
{
  file=$1
  scalar_name=$2
  scalar_val=$3
  sed -i "s|^\s*${scalar_name}\s*:\s*.*$|${scalar_name}: ${scalar_val}|" $file
}

runSims() {
  name=$1
  mechFile=$2
  thermFile=$3
  sparseThresh=$4
  ignitionDelayTime=$5
  dodense=$6

  ignitionDelayTime=`echo ${ignitionDelayTime} | sed 's/[eE]/\\*10\\^/' | sed 's/+//'`
  simTime=$(echo "scale=20; ${ignitionDelayTime} * 2" | bc)


  initFuelMassFracs=`grep -v "^#" inputs/${name}_fracs.log | egrep -iv '^\<O2\>|^\<N2\>' | cut -f1-2 | sed 's/$/,/' | sed 's/\s\s*/: /g'`
  initFuelMassFracs=`echo $initFuelMassFracs`  #strip newlines
  initOxidMassFracs=`grep -v "^#" inputs/${name}_fracs.log | egrep -i  '^\<O2\>|^\<N2\>' | cut -f1-2 | sed 's/$/,/' | sed 's/\s\s*/: /g'`
  initOxidMassFracs=`echo $initOxidMassFracs`  #strip newlines

  name=${name}_test

  #infile=run.yml
  infile=`mktemp --tmpdir=./`
  cp inputs/base.yml $infile

  setYMLScalar $infile mechFile inputs/$mechFile
  setYMLScalar $infile thermFile inputs/$thermFile
  setYMLScalar $infile reactorTime $simTime
  setYMLScalar $infile precThresh $sparseThresh
  setYMLScalar $infile fuelComp "{ $initFuelMassFracs }"
  setYMLScalar $infile oxidizerComp "{ $initOxidMassFracs }"
  setYMLScalar $infile mechLogFile outputs/${name}.cklog
  setYMLScalar $infile outFile outputs/${name}_mr_run.log
  
  ork_exe=../zerork_reactor.x


  #CPU Dense 
  if [ "x"$dodense == "xY" ]
  then
    echo "Running ${name} CPU Dense"
    setYMLScalar $infile cuda off
    setYMLScalar $infile linearSolver DirectDenseDVD
    setYMLScalar $infile nReactors 32 
    setYMLScalar $infile nMatrixReactors 32 
    setYMLScalar $infile outFile outputs/${name}_mr_run_cpu_dense.log
    #mpirun -np 1 $ork_exe $infile > outputs/log_${name}_cpu_dense
    $ork_exe $infile > outputs/log_${name}_cpu_dense
  fi

  #CPU Sparse
  echo "Running ${name} CPU Sparse"
  setYMLScalar $infile cuda off
  setYMLScalar $infile linearSolver IterativeSparse 
  setYMLScalar $infile nReactors 32 
  setYMLScalar $infile nMatrixReactors 32 
  setYMLScalar $infile outFile outputs/${name}_mr_run_cpu_sparse.log
  #mpirun -np 1 $ork_exe $infile > outputs/log_${name}_cpu_sparse
  $ork_exe $infile > outputs/log_${name}_cpu_sparse
  rm $infile
}


nm=h2
mf=h2_v1b_mech.txt
tf=h2_v1a_therm.txt
st=3.20e-5
idt=8.2565998e-07

runSims $nm $mf $tf $st $idt Y

nm=dme
mf=dme_24_mech.txt
tf=dme_24_therm.txt
st=6.40e-5
idt=1.1e-05

runSims $nm $mf $tf $st $idt Y

nm=nc7h16-skel
mf=heptanesymp159_mec.txt
tf=heptanesymp_therm.txt
st=1.28e-4
idt=1.35e-05

runSims $nm $mf $tf $st $idt N



