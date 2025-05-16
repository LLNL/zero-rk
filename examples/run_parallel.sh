
examples_dir=$PWD

for d in \
  perturbAFactorGSA \
  variable_volume_gsa \
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
  surrogate_optimizer/idt_sweep
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

