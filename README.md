#  Zero-RK : Zero Order Reaction Kinetics

Zero-RK is a software package that simulates chemically reacting systems in a computationally efficient manner. The fundamental advance embodied in Zero-RK is the numerical approach, which results in orders-of-magnitude reduction in simulation time while maintaining the accuracy of the results.

Authors
----------------
Zero-RK was created by Matthew McNenly, mcnenly1@llnl.gov, with significant contributions from Russell Whitesides whitesides1@llnl.gov.

Build Notes
----------------

To build all the dependencies and the two Zero-RK example applications execute
the script build_all.sh in this directory. If there are any build errors like

/usr/bin/ld: cannot find -llapack
/usr/bin/ld: cannot find -lblas

this indicates that the makefile commands do not have the correct paths to 
certain common system libraries.  The file build.config in this directory 
may need to be updated to correct the path information for the following
system libraries:

1. LAPACK (e.g. liblapack.a or liblapack.so)
2. BLAS (e.g. liblas.a or liblas.so)

Applications
----------------

In the current release, two applications are included:

1. [Constant Volume Well-Stirred Reactor](https://github.com/llnl/zero-rk/blob/master/zerork/applications/constVolumeWSR)

2. [Zero-RK Reactor Interface Library](https://github.com/llnl/zero-rk/blob/master/zerork_reactor)

License
----------------

Zero-RK is distributed under the terms of the BSD-3-Clause license.

All new contributions must be made under the same license.

See LICENSE and NOTICE for details.

SPDX-License-Identifier: (BSD-3-Clause)

LLNL-CODE-779961
