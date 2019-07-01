# Constant Volume Well-stirred Reactor

Three example test cases are included to compute constant volume
well-stirred reactor ignition delay times.  These cases compute the
ignition delay time for a range of initial temperatures for three mechanisms.

The test case sub-directories each contain an input file that is used as
follows (the hydrogen mechanism is shown below):

```sh
cd examples/hydrogen
../../constVolumeWSR.x hydrogen_TSweep.inp
```

This will create three output files:

| filename | file contents |
| --- | --- |
| `h2_v1b_mech.log`       | mechanism parser log file |
| `hydrogen_TSweep.dat`   | ignition delay times |
| `hydrogen_TSweep.thist` | time history (print frequency currently set to 1 sec in input file) |

Reference examples of these files are also included with `*.gold` suffix for
comparison and verification.

There is another constant volume well-stirred reactor executable that uses
the traditional dense ODE solve in cvode with lapack, which can be used as
an additional check. It uses Zero-RK for the computation of the time
derivatives, but relies on the automatic divided difference computation
of the dense Jacobians. It uses the same input files as found in the 
examples sub-directories.  For large mechanisms, the performance is
significantly slower than Zero-RK's adaptive preconditioner approach. 
Specifically, to test performance with this application, run the
following in the `examples/iso_octane` directory:

```sh
    ../../direct_CVDense/cv_wsr_direct_CVDense.x iso_octane_TSweep.inp
```
    
Reference examples using the `hydrogen` and `iso_octane` mechanisms for 
dense simulations are included with `*.gold.dense` suffix for comparison and 
verification.

