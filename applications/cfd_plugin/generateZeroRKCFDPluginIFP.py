#!/usr/bin/env python

import os,sys

spify_parser_name = "ZeroRKCFDPluginIFP"

spify_parser_params = []

spify_parser_params.append(
{
    'name':'mechanism_parsing_log',
    'type':'string',
    'shortDesc' : "Mechanism Parser Log File",
    'defaultValue' : '/dev/null'
}
)

spify_parser_params.append(
{
    'name':'reactor_timing_log',
    'type':'string',
    'shortDesc' : "Log file for reactor solutions timing statistics",
    'defaultValue' : '/dev/null'
}
)

spify_parser_params.append(
{
    'name':'n_reactors_max',
    'type':'int',
    'shortDesc' : "Maximum number of reactors to solve",
    'boundMin': 1,
    'boundMax': 8192,
    'defaultValue': 8192
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
    'defaultValue':1000000
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
    'discreteValues': [0,1]
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
    'name':"load_balance",
    'type':'int',
    'shortDesc' : "Do internal load balance",
    'defaultValue' : 1,
    'discreteValues': [0,1]
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
    'defaultValue' : 400.0
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
    'name':"dump_reactors",
    'type':'int',
    'shortDesc' : "Dump all reactors to files",
    'defaultValue' : 0,
    'discreteValues': [0,1]
}
)


from SpifyParserGenerator import SpifyParserGenerator as spg

spg().generate(spify_parser_name,spify_parser_params)
#spg().make_master_file(spify_parser_name,spify_parser_params)


