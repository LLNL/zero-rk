#!/usr/bin/env python

import os

from SpifyParserGenerator import SpifyParserGenerator as spg

#Name your parser
spify_parser_name = "ThermoCheckIFP"

spify_parser_params = []

#Specify parameters
#ork_cvIDT_multireactor parms

spify_parser_params.append(
{
    'name':'mech_file',
    'type':'string',
    'shortDesc' : 'Chemkin format mechansim file',
    'defaultValue' : os.devnull
}
)

spify_parser_params.append(
{
    'name':'therm_file',
    'type':'string',
    'shortDesc' : 'Chemkin format thermodynamics file'
}
)

spify_parser_params.append(
{
    'name':'repair_file',
    'type':'string',
    'shortDesc' : 'Repaired Chemkin format thermodynamics file'
}
)

spify_parser_params.append(
{
    'name' : 'add_comments',
    'type' : 'bool',
    'longDesc' : 'If set to true [1], then the program will include the comments about the maximum difference in Cp/R, H/RT and S/R for each species in [repair_file]',
    'defaultValue' : 1
}
)

spify_parser_params.append(
{
    'name' : 'change_file',
    'type' : 'string',
    'longDesc' : 'File with tables suitable for plotting the properties Cp/R, H/RT and S/R to inspect the repair changes',
    'defaultValue' : os.devnull
}
)

spify_parser_params.append(
{
    'name':'detailed_table_file',
    'type':'string',
    'shortDesc' : 'File with Detailed Table of Thermodyanic Info',
    'defaultValue' : ""
}
)

spify_parser_params.append(
{
    'name':'isomer_plot_file',
    'type':'string',
    'shortDesc' : 'File with Detailed Table of Thermodyanic Info',
    'defaultValue' : ""
}
)

spify_parser_params.append(
{
    'name' : 'find_best_Tmatch',
    'type' : 'bool',
    'longDesc' : 'If set to true [1], then the program will find try to find a new matching temperature where the Cp/R jump is minimized. If set to false [0], then the program will keep the original matching temperature for the species',
    'defaultValue' : 0
}
)

spify_parser_params.append(
{
    'name' : 'fix_Tmatch',
    'type' : 'bool',
    'longDesc' : 'If set to true [1], then the program will use the global matching temperature specified on the line following the THERMO keyword for all species',
    'defaultValue' : 0
}
)

spify_parser_params.append(
{
    'name' : 'use_min_species',
    'type' : 'bool',
    'longDesc' : 'If set to true [1], then the program will only repair and write the first appearance of species in the thermodynamics file [therm_file], and the subset of species in the mechanism file if provided [mech_file]',
    'defaultValue' : 0
}
)

spify_parser_params.append(
{
    'name':'repair_atol',
    'type':'double',
    'longDesc' : 'Absolute tolerance of the change in the repaired Cp/R, H/RT and S/R beyond which the thermodynamics differences are reported in the [change_file]',
    'defaultValue' : 0.1
}
)

spify_parser_params.append(
{
    'name':'retain_temperature',
    'type':'double',
    'longDesc' : 'Temperature at which the properties from the original thermodynamics are retained, for example to keep the same standard state properties',
    'defaultValue' : 273.15
}
)

spify_parser_params.append(
{
    'name':'num_repair_points',
    'type':'int',
    'longDesc' : 'Number of points to minimize the difference in the repaired thermodynamic properties using the constrained linear least squares method',
    'defaultValue' : 100
}
)

spify_parser_params.append(
{
    'name':'num_test_points',
    'type':'int',
    'longDesc' : 'Number of points to test when checking the difference in the repaired thermodynamic properties',
    'defaultValue' : 1000
}
)

spify_parser_params.append(
{
    'name':'num_print_points',
    'type':'int',
    'longDesc' : 'Number of points printed per species in the tables recorded in the [change_file]',
    'defaultValue' : 200
}
)






#Generate parser code
spg().generate(spify_parser_name,spify_parser_params)


#Done.
