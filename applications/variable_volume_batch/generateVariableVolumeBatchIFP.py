#!/usr/bin/env python

import os

from SpifyParserGenerator import SpifyParserGenerator as spg

#Name your parser
spify_parser_name = "VariableVolumeBatchIFP"

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
    'name':'initialConditionsFile',
    'type':'string',
    'longDesc' : "File containing initial conditions. (Column 1) Fuel mole fraction; (2) Ethanol mole fraction; (3) O2 mole fraction; (4) N2 mole fraction; (5) AR mole fraction; (6) Initial temperature [C]; (7) Inital pressure [bar]; and (8) volume history file [csv format t, vol]"
}
)

spify_parser_params.append(
{
    'name':'volumeFilePath',
    'type':'string',
    'shortDesc' : "Path to directory containing volume histories",
    'defaultValue' : './'
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
    'name':'timeHistoryOutputFileBase',
    'type':'string',
    'shortDesc' : "Base name for time history output files",
    'defaultValue' :  ""
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
    'name':'outputFile',
    'type':'string',
    'shortDesc' : "Output summary file"
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
    'name':'oxidizerSpecies',
    'type':'v_string',
    'shortDesc' : "Species in oxidizer to be read from initial conditions file"
}
)

spify_parser_params.append(
{
    'name':'fuelComp',
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
