# Zero-RK CFD Plugin Library

The `zerork_cfd_plugin_tester` application implements a simple interface for 
calculation of chemical source terms that can be coupled to computational
fluid dynamics simulations.

An [test case] is provided to verify that the code is working correctly.

The test includes a script to exercise the `zerork_reactor` application
which simultaneously solves multiple constant volume well-stirred reactors.

The results of the test case simulations are stored in a subdirectory named
`outputs`.  These results can then be compared against the reference data
in the `outputs.gold` subdirectory for verification.

[test case]: https://github.com/llnl/zero-rk/blob/master/applications/cfd_plugin/test

