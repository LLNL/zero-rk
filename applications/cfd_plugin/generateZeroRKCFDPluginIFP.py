#!/usr/bin/env python

import os,sys

spify_parser_name = "ZeroRKCFDPluginIFP"

spify_parser_params = []

spify_parser_params.append(
{
    'name':'mechanism_parsing_log',
    'type':'string',
    'shortDesc' : "Mechanism Parser Log File",
    'defaultValue' : os.devnull
}
)

spify_parser_params.append(
{
    'name':'reactor_timing_log',
    'type':'string',
    'shortDesc' : "Log file for reactor solutions timing statistics",
    'defaultValue' : os.devnull
}
)

spify_parser_params.append(
{
    'name':'n_reactors_max',
    'type':'int',
    'shortDesc' : "Maximum number of reactors to solve",
    'boundMin': 1,
    'defaultValue': 2048
}
)

spify_parser_params.append(
{
    'name':'n_reactors_min',
    'type':'int',
    'shortDesc' : "Minimum number of reactors to solve",
    'boundMin': 1,
    'defaultValue': 256
}
)

spify_parser_params.append(
{
    'name':'max_dt',
    'type':'double',
    'shortDesc' : "Maximum Internal Integrator Step Size",
    'defaultValue' : 0.05
}
)

spify_parser_params.append(
{
    'name' :'max_steps',
    'type':'int',
    'shortDesc' : "Maximum Number of Integrator Steps",
    'defaultValue': 5000
}
)

spify_parser_params.append(
{
    'name':'relative_tolerance',
    'type':'double',
    'shortDesc':"Relative Integrator Tolerance",
    'defaultValue':1.0e-8
}
)

spify_parser_params.append(
{
    'name':"absolute_tolerance",
    'type':'double',
    'shortDesc':"Absolute Integrator Tolerance",
    'defaultValue':1.0e-20
}
)

spify_parser_params.append(
{
    'name':"abstol_dens",
    'type':'int',
    'shortDesc':"Use density based absolute tolerance",
    'defaultValue': 0
}
)


spify_parser_params.append(
{
    'name':"reference_temperature",
    'type':'double',
    'shortDesc' : "Ignition Delat Metric: Maximum species concentration",
    'defaultValue' : 1.0
}
)

spify_parser_params.append(
{
    'name':"preconditioner_threshold",
    'type':'double',
    'shortDesc' : "Preconditioner threshold values"
}
)

spify_parser_params.append(
{
    'name':"eps_lin",
    'type':'double',
    'shortDesc' : "??",
    'defaultValue' : 0.1
}
)

spify_parser_params.append(
{
    'name':"nonlinear_convergence_coeff",
    'type':'double',
    'shortDesc' : "??",
    'defaultValue' : 0.05
}
)

spify_parser_params.append(
{
    'name':"integrator",
    'type':'int',
    'shortDesc' : "Integrator",
    'defaultValue' : 0,
    'discreteValues': [0,1,2,3]
}
)

spify_parser_params.append(
{
    'name':"dense",
    'type':'int',
    'shortDesc' : "Dense Matrix",
    'defaultValue' : 1,
    'discreteValues': [0,1]
}
)

spify_parser_params.append(
{
    'name':"analytic",
    'type':'int',
    'shortDesc' : "Analytic Jacobian",
    'defaultValue' : 1,
    'discreteValues': [0,1]
}
)

spify_parser_params.append(
{
    'name':"iterative",
    'type':'int',
    'shortDesc' : "Iterative (1) or Direct (0)",
    'defaultValue' : 1,
    'discreteValues': [0,1]
}
)

spify_parser_params.append(
{
    'name':"gpu",
    'type':'int',
    'shortDesc' : "GPU Option",
    'defaultValue' : 0,
    'discreteValues': [0,1]
}
)

spify_parser_params.append(
{
    'name':"initial_gpu_multiplier",
    'type':'double',
    'shortDesc' : "Initial value for how much more work to give GPU",
    'defaultValue' : 20.0,
    'boundMin': 1.0
}
)

spify_parser_params.append(
{
    'name':"load_balance",
    'type':'int',
    'shortDesc' : "Which load balance method.",
    'defaultValue' : 1,
    'discreteValues': [0,1,2]
}
)

spify_parser_params.append(
{
    'name':"load_balance_noise",
    'type':'int',
    'shortDesc' : "Parameter to improve load balance",
    'defaultValue' : 0,
    'boundMin': 0
}
)

spify_parser_params.append(
{
    'name':"load_balance_mem",
    'type':'int',
    'shortDesc' : "Which load-balancing memory option to use. 0 uses less memory but is slower than 1",
    'defaultValue' : 1,
    'boundMin': 0
}
)

spify_parser_params.append(
{
    'name':"reactor_weight_mult",
    'type':'int',
    'shortDesc' : "Another parameter to improve load balance",
    'defaultValue' : 1,
    'boundMin': 1
}
)

spify_parser_params.append(
{
    'name':"sort_reactors",
    'type':'int',
    'shortDesc' : "Sort reactors by cost",
    'defaultValue' : 1,
    'discreteValues': [0,1]
}
)

spify_parser_params.append(
{
    'name':"verbosity",
    'type':'int',
    'shortDesc' : "Output verbosity",
    'defaultValue' : 4,
    'discreteValues': [0,1,2,3,4]
}
)

spify_parser_params.append(
{
    'name':"constant_volume",
    'type':'int',
    'shortDesc' : "Constant volume or constant pressure",
    'defaultValue' : 1,
    'discreteValues': [0,1]
}
)

spify_parser_params.append(
{
    'name':"delta_temperature_ignition",
    'type':'double',
    'shortDesc' : "",
    'defaultValue' : 0.0,
    'boundMin': 0.0,
    'boundMax': 100000.0
}
)

spify_parser_params.append(
{
    'name':"min_mass_fraction",
    'type':'double',
    'shortDesc' : "",
    'defaultValue' : 1.0e-30
}
)

spify_parser_params.append(
{
    'name':"step_limiter",
    'type':'double',
    'shortDesc' : "",
    'defaultValue' : 1.0e22
}
)

spify_parser_params.append(
{
    'name':"always_solve_temperature",
    'type':'int',
    'shortDesc' : "Whether to enable/disable constant temperature reactor type",
    'defaultValue' : 1,
    'discreteValues': [0,1]
}
)

spify_parser_params.append(
{
    'name':"solve_temperature_threshold",
    'type':'double',
    'shortDesc' : "Temperature delta from previous step which disables constant-temperature solve (if enabled).",
    'defaultValue' : 2.0,
    'boundMin': 0
}
)

spify_parser_params.append(
{
    'name':"dump_reactors",
    'type':'int',
    'shortDesc' : "Dump all reactors to files",
    'defaultValue' : 0,
    'discreteValues': [0,1]
}
)

spify_parser_params.append(
{
    'name':"dump_failed_reactors",
    'type':'int',
    'shortDesc' : "Dump failed reactors to files",
    'defaultValue' : 0,
    'discreteValues': [0,1]
}
)

spify_parser_params.append(
{
    'name':"output_performance_log",
    'type':'int',
    'shortDesc' : "Write performance log",
    'defaultValue' : 1,
    'discreteValues': [0,1]
}
)

spify_parser_params.append(
{
    'name':"cvode_num_retries",
    'type':'int',
    'shortDesc' : "Number of times to retry after CVODE failure",
    'defaultValue' : 3,
    'boundMin':  0,
    'boundMax':  20
}
)

spify_parser_params.append(
{
    'name':"cvode_retry_absolute_tolerance_adjustment",
    'type':'double',
    'shortDesc' : "Mulitplicative factor to change absolute tolerance on retrying after CVODE failure",
    'defaultValue' : 0.1,
    'boundMin':  1.0e-6,
    'boundMax':  1.0
}
)

spify_parser_params.append(
{
    'name':"cvode_retry_relative_tolerance_adjustment",
    'type':'double',
    'shortDesc' : "Mulitplicative factor to change relative tolerance on retrying after CVODE failure",
    'defaultValue' : 1,
    'boundMin':  1.0e-6,
    'boundMax':  1.0
}
)


from SpifyParserGenerator import SpifyParserGenerator as spg

spg().generate(spify_parser_name,spify_parser_params)
#spg().make_master_file(spify_parser_name,spify_parser_params)


