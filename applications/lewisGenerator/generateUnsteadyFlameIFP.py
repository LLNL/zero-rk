#!/usr/bin/env python

import os

from SpifyParserGenerator import SpifyParserGenerator as spg

#Name your parser
spify_parser_name = "UnsteadyFlameIFP"

spify_parser_params = []

#Specify parameters
spify_parser_params.append(
{
    'name':'mech_file',
    'type':'string',
    'shortDesc' : "Chemkin Format Mechansim File [input]"
}
)

spify_parser_params.append(
{
    'name':'therm_file',
    'type':'string',
    'shortDesc' : "Chemkin Format Thermodynamics File [input]"
}
)

spify_parser_params.append(
{
    'name':'trans_file',
    'type':'string',
    'shortDesc' : "Transport File [input]"
}
)

spify_parser_params.append(
{
    'name':'transport_model',
    'type':'string',
    'shortDesc' : "Transport model [input]"
}
)

spify_parser_params.append(
{
    'name':'log_file',
    'type':'string',
    'shortDesc' : "Mechanism Parser Log File [output]",
    'defaultValue' : os.devnull
}
)


spify_parser_params.append(
{
    'name':'pressure',
    'type':'double',
    'shortDesc' : "Pressure of the flame system [Pa]",
    'defaultValue' : 1.01325e5
}
)

spify_parser_params.append(
{
    'name':'inlet_fuel_comp',
    'type':'m_string_double',
    'shortDesc' : "Fuel composition mole fraction map for inlet"
}
)

spify_parser_params.append(
{
    'name':'inlet_oxidizer_comp',
    'type':'m_string_double',
    'shortDesc' : "Oxidizer composition mole fraction map for inlet"
}
)

spify_parser_params.append(
{
    'name':'inlet_full_comp',
    'type':'m_string_double',
    'shortDesc' : "Full composition mole fraction map for inlet"
}
)

spify_parser_params.append(
{
    'name':"ref_temperature",
    'type':'double',
    'shortDesc' : "Reference temperature to normalize the ODE dT/dt equation",
    'defaultValue' : 1000.0
}
)

spify_parser_params.append(
{
    'name':"inlet_temperature",
    'type':'double',
    'shortDesc' : "Inlet temperature [K]",
    'defaultValue' : 1800.0
}
)

spify_parser_params.append(
{
    'name':"inlet_phi",
    'type':'double',
    'shortDesc' : "Inlet equivalence ratio (F/A)/(F/A)_{st} [-]",
    'defaultValue': 1.0
}
)

spify_parser_params.append(
{
    'name':"egr",
    'type':'double',
    'shortDesc' : "Inlet EGR value [-]",
    'defaultValue': 0.0
}
)


spify_parser_params.append(
{
    'name':"use_equilibrium",
    'type':'bool',
    'shortDesc' : "Flag that when set to true [y] computes the H-P equilibrium composition to evaluate Lewis number. Inlet values used otherwise",
    'defaultValue': 0
}
)



#Generate parser code
spg().generate(spify_parser_name,spify_parser_params)


#Done.
