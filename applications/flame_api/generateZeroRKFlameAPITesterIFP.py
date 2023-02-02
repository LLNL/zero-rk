#!/usr/bin/env python

import os

from SpifyParserGenerator import SpifyParserGenerator as spg

#Name your parser
spify_parser_name = "ZeroRKFlameAPITesterIFP"

spify_parser_params = []

#Specify parameters
spify_parser_params.append(
{
    'name':'mechanism_file',
    'type':'string',
    'shortDesc' : "Chemkin Format Mechansim File"
}
)

spify_parser_params.append(
{
    'name':'therm_file',
    'type':'string',
    'shortDesc' : "Chemkin Format Thermodynamics File"
}
)

spify_parser_params.append(
{
    'name':'transport_file',
    'type':'string',
    'shortDesc' : "Transport File"
}
)

spify_parser_params.append(
{
    'name':'transport_model',
    'type':'string',
    'shortDesc' : "Transport model",
    'discreteValues' : ["ConstantLewis", "MixAvg", "MixAvgSoret"]
}
)


spify_parser_params.append(
{
    'name':'flame_start_profile',
    'type':'string',
    'longDesc' : "ASCII Input File containing grid, temperature and mass fractions",
}
)

spify_parser_params.append(
{
    'name':'flame_end_profile',
    'type':'string',
    'longDesc' : "Output CSV file containing grid, temperature and mass fractions",
}
)

spify_parser_params.append(
{
    'name':'pressure',
    'type':'double',
    'shortDesc' : "Pressure of the flame system [Pa]"
}
)

spify_parser_params.append(
{
    'name':'flame_speed',
    'type':'double',
    'shortDesc' : "Flame speed of initial solution"
}
)

spify_parser_params.append(
{
    'name':'relative_tolerance',
    'type':'double',
    'shortDesc':"Relative Integrator Tolerance",
    'defaultValue': 1.0e-2
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
    'name':"verbosity",
    'type':'int',
    'shortDesc':"Verbosity level",
    'defaultValue': 0
}
)

spify_parser_params.append(
{
    'name':"pseudo_unsteady",
    'type':'int',
    'shortDesc':"Run psuedo_unsteady at first",
    'defaultValue': 0
}
)

spify_parser_params.append(
{
    'name':"pseudo_unsteady_dt",
    'type':'double',
    'shortDesc' : "Initial pseudo-unsteady time step",
    'defaultValue': 1.0e-3
}
)

spify_parser_params.append(
{
    'name':"pseudo_unsteady_max_iterations",
    'type':'int',
    'shortDesc':"Number of iterations for psuedo unsteady solution",
    'defaultValue': 80
}
)

spify_parser_params.append(
{
    'name':"pseudo_unsteady_time",
    'type':'double',
    'shortDesc' : "Time duration for pseudo-unsteady solution",
    'defaultValue': 0.05
}
)

spify_parser_params.append(
{
    'name':"step_limiter",
    'type':'double',
    'shortdesc' : "reaction rate limiter",
    'defaultvalue': 1.0e+300
}
)

spify_parser_params.append(
{
    'name':'integrator_type',
    'type':'int',
    'longDesc' : "Integrator type: (0, 1 = unsupported, 2 = Exact Jacobian with SuperLU solver, 3 = Approximate Jacobian with SuperLU + LAPACK solvers)",
    'defaultValue' : 3
}
)

spify_parser_params.append(
{
    'name':'convective_scheme_type',
    'type':'int',
    'longDesc' : "Convective scheme type: (0 = First order upwind, 1 = Second order upwind, 2 = Second order centered)",
    'defaultValue' : 2
}
)

spify_parser_params.append(
{
    'name':"reference_temperature",
    'type':'double',
    'shortDesc' : "Reference temperature to normalize the ODE dT/dt equation",
    'defaultValue' : 1000.0
}
)

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
    'name':'transport_parsing_log',
    'type':'string',
    'shortDesc' : "Transport Parser Log File",
    'defaultValue' : os.devnull
}
)

spify_parser_params.append(
{
    'name':'store_jacobian',
    'type':'int',
    'longDesc' : "Flag that when set non-zero stores the jacobian data for each grid point when computed by user (integrator_type == 3)",
    'defaultValue' : 1
}
)

spify_parser_params.append(
{
    'name':'steady_max_iterations',
    'type':'int',
    'shortDesc' : "Maximum number of KINsol iterations during steady phase",
    'defaultValue' : 100
}
)

spify_parser_params.append(
{
    'name':'max_subiterations',
    'type':'int',
    'shortDesc' : "Maximum number of iterations between jacobian updates",
    'defaultValue' : 10
}
)

#Generate parser code
spg().generate(spify_parser_name,spify_parser_params)


#Done.
