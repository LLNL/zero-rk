#!/usr/bin/env python

import os

from SpifyParserGenerator import SpifyParserGenerator as spg

#Name your parser
spify_parser_name = "VariableVolumeIFP"

spify_parser_params = []

#Specify parameters
#ork_cvIDT_multireactor parms

spify_parser_params.append(
{
    'name':'mechFile',
    'type':'string',
    'shortDesc' : "Chemkin Format Mechansim File"
}
)

spify_parser_params.append(
{
    'name':'thermFile',
    'type':'string',
    'shortDesc' : "Chemkin Format Thermodynamics File"
}
)

spify_parser_params.append(
{
    'name':'mechLogFile',
    'type':'string',
    'shortDesc' : "Mechanism Parser Log File",
    'defaultValue' : os.devnull
}
)

spify_parser_params.append(
{
    'name':'volumeFile',
    'type':'string',
    'longDesc' : "File containing the volume time history (column 1 = time [s], column 2 = volume [m^3])",
}
)

spify_parser_params.append(
{
    'name':'strokeMult',
    'type':'double',
    'longDesc' : "Stroke multiplier [-] to rescale the volume time history (preserves the minimum volume)",
    'defaultValue' : 1.0
}
)

spify_parser_params.append(
{
    'name':'volumeMult',
    'type':'double',
    'shortDesc' : "Volume multiplier [-] to rescale the volume time history",
    'defaultValue' : 1.0
}
)

spify_parser_params.append(
{
    'name':'thistFile',
    'type':'string',
    'shortDesc' : "Time history output file"
}
)

spify_parser_params.append(
{
    'name':'timeHistoryOutputStepPeriod',
    'type':'int',
    'shortDesc' : "Time history output printing frequency in steps",
    'defaultValue' :  "10"
}
)

spify_parser_params.append(
{
    'name':'timeHistoryOutputMinimumTimeStep',
    'type':'double',
    'shortDesc' : "Time history output printing frequency in simulation time",
    'defaultValue' :  "1.0e-5"
}
)

spify_parser_params.append(
{
    'name':'timeHistoryOutputHighPrecision',
    'type':'bool',
    'shortDesc' : "Use high precision output for time history",
    'defaultValue' :  1
}
)


spify_parser_params.append(
{
    'name':'initTime',
    'type':'double',
    'shortDesc' : "Initial ODE system time"
}
)

spify_parser_params.append(
{
    'name':'finalTime',
    'type':'double',
    'shortDesc' : "Final ODE system time"
}
)

spify_parser_params.append(
{
    'name':'maxDtInternal',
    'type':'double',
    'shortDesc' : "Maximum Internal Integrator Step Size",
    'defaultValue' : 0.01
}
)

spify_parser_params.append(
{
    'name' :'maxSteps',
    'type':'int',
    'shortDesc' : "Maximum Number of Integrator Steps",
    'defaultValue':1000000
}
)

spify_parser_params.append(
{
    'name':'relTol',
    'type':'double',
    'shortDesc':"Relative Integrator Tolerance",
    'defaultValue':1.0e-8
}
)

spify_parser_params.append(
{
    'name':"absTol",
    'type':'double',
    'shortDesc':"Absolute Integrator Tolerance",
    'defaultValue':1.0e-20
}
)

spify_parser_params.append(
{
    'name':"refTemp",
    'type':'double',
    'shortDesc' : "Reference temperature to normalize the ODE state variables",
    'defaultValue' : 1000.0
}
)

spify_parser_params.append(
{
    'name':"refMoles",
    'type':'double',
    'shortDesc' : "Reference moles to normalize the ODE state variables",
    'defaultValue' : 1.0
}
)

spify_parser_params.append(
{
    'name':'initComp',
    'type':'m_string_double',
    'shortDesc' : "Initial composition map"
}
)

spify_parser_params.append(
{
    'name':'compType',
    'type':'string',
    'shortDesc' : "Composition type (MOLE_FRACTION or MASS_FRACTION)",
    'defaultValue' : "MOLE_FRACTION"
}
)

spify_parser_params.append(
{
    'name':"initTemp",
    'type':'double',
    'shortDesc' : "Initial temperature [K]",
}
)

spify_parser_params.append(
{
    'name':"initPres",
    'type':'double',
    'shortDesc' : "Initial Pressure [Pa]",
}
)

spify_parser_params.append(
{
    'name':"maxOrder",
    'type':'int',
    'shortDesc' : "Maximum BDF",
    'defaultValue' : 5
}
)

spify_parser_params.append(
{
    'name':"maxNumNonLinIters",
    'type':'int',
    'shortDesc' : "Maximum number of nonlinear interations",
    'defaultValue' : 3
}
)

spify_parser_params.append(
{
    'name':"cvEpsLin",
    'type':'double',
    'shortDesc' : "Linear solver tolerance multiplier",
    'defaultValue' : 0.1
}
)

spify_parser_params.append(
{
    'name':"cvNlConvCoeff",
    'type':'double',
    'shortDesc' : "Non-linear solver tolerance multiplier",
    'defaultValue' : 0.05
}
)

spify_parser_params.append(
{
    'name':"Jacobian",
    'type':'string',
    'shortDesc' : "Type of Jacobian solution method",
    'defaultValue' : "SPARSE_ANALYTIC",
    'discreteValues': ["SPARSE_ANALYTIC","DENSE_ANALYTIC","DENSE_NUMERIC"]
}
)

spify_parser_params.append(
{
    'name':"negativeConcentrationMitigation",
    'type':'int',
    'shortDesc' : "Prevent explosive growth of negative mass fractions",
    'defaultValue' : 0
}
)

spify_parser_params.append(
{
    'name':"volumeInterpolationMethod",
    'type':'string',
    'shortDesc' : "Interpolation method for volume time history",
    'defaultValue' : "LINEAR_CLIPPED",
    'discreteValues': ["LINEAR_CLIPPED","CUBIC_CLIPPED"]
}
)

#Generate parser code
spg().generate(spify_parser_name,spify_parser_params)


#Done.
