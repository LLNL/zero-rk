
examples_dir=$PWD

for d in \
  constVolumeWSR \
  constVolumeWSR_TLA \
  constVolumePSR \
  idt_diagnostic \
  perturbAFactor \
  thermo_check \
  variable_volume/cv \
  variable_volume/rcm \
  variable_volume_batch \
  cfd_plugin_tester \
  rate_optimization
do
  cd $d
  echo "Running test $d ..."
  sh run.sh >& run.log
  if [ ! "$?" -eq 0 ]
  then
    echo ""
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "   Example: $d failed.     "
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo ""
  fi
  cd $examples_dir
done

