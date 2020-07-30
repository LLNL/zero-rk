
examples_dir=$PWD

for d in \
  constVolumeWSR/hydrogen \
  constVolumeWSR/hydrogen_yml \
  constVolumeWSR/iso_octane \
  diffusion_steady_flame_solver \
  diffusion_unsteady_flame_solver \
  idt_diagnostic \
  lewisGenerator \
  flame_speed/steady \
  flame_speed/unsteady \
  perturbAFactor \
  perturbAFactorGSA \
  premixed_steady_flame_solver \
  premixed_unsteady_flame_solver \
  thermo_check \
  variable_volume/cv \
  variable_volume/rcm \
  variable_volume_batch \
  variable_volume_gsa \
  cfd_plugin_tester
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

