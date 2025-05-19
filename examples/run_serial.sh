
examples_dir=$PWD

for d in \
  constVolumeWSR \
  constVolumeWSR_TLA \
  constVolumePSR \
  idt_diagnostic \
  thermo_check \
  perturbAFactor \
  variable_volume/cv \
  variable_volume/rcm \
  variable_volume_batch \
  lewisGenerator \
  cfd_plugin_tester \
  rate_optimization \
  surrogate_optimizer/diesel \
  flame_api_tester \
  #surrogate_optimizer/tprf \
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

