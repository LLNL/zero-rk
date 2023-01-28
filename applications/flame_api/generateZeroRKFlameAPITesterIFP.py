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
    'defaultValue':1.0e-2
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
    'name':"step_limiter",
    'type':'double',
    'shortDesc' : "Reaction rate limiter",
    'defaultValue': 1.0e+300
}
)


#Generate parser code
spg().generate(spify_parser_name,spify_parser_params)


#Done.
