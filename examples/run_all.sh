
examples_dir=$PWD

for d in \
  constVolumeWSR \
  constVolumeWSR_TLA \
  constVolumePSR \
  idt_diagnostic \
  perturbAFactor \
  perturbAFactorGSA \
  variable_volume/cv \
  variable_volume/rcm \
  variable_volume_batch \
  variable_volume_gsa \
  lewisGenerator \
  premixed_steady_flame_solver \
  premixed_unsteady_flame_solver \
  diffusion_steady_flame_solver \
  diffusion_unsteady_flame_solver \
  counterflow_steady_flame_solver \
  counterflow_unsteady_flame_solver \
  flame_speed/steady \
  flame_speed/unsteady \
  counterflow_flame_sweeps/diffusion_steady \
  counterflow_flame_sweeps/diffusion_unsteady \
  counterflow_flame_sweeps/premixed_steady \
  counterflow_flame_sweeps/premixed_unsteady \
  thermo_check \
  cfd_plugin_tester \
  rate_optimization \
  flame_api_tester
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

