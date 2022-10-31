#!/usr/bin/env python

import os

from SpifyParserGenerator import SpifyParserGenerator as spg

#Name your parser
spify_parser_name = "VariableVolumeGSAIFP"

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
    'name':'gsaMatrixFile',
    'type':'string',
    'shortDesc' : "Input file containing the multiplicative A-Factor perturbations"
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
    'name':'refHeatRelease',
    'type':'double',
    'longDesc' : "Reference heat release [J] to compute the fraction of chemical heat release (negative is exothermic)",
    'defaultValue' : -1.0
}
)

spify_parser_params.append(
{
    'name':'useRefHeatRelease',
    'type':'bool',
    'longDesc' : "Use the prescribed refHeatRelease value (1=true, 0=false), othherwise use the heat release computed for the unperturbed simulation",
    'defaultValue' : 0
}
)

spify_parser_params.append(
{
    'name':"burnFraction",
    'type':'v_double',
    'shortDesc' : "Fraction(s) of chemical heat release at which the time is reported",
    'defaultValue' : [ 0.1, 0.5, 0.9 ]
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
    'name':'thistMinDt',
    'type':'double',
    'shortDesc' : "Minimum time change to print state to time history file",
    'defaultValue' : 1.0e-6
}
)

spify_parser_params.append(
{
    'name':'thistStepPeriod',
    'type':'int',
    'longDesc' : "Number of internal time steps before printing the state to the time history file ",
    'defaultValue' : 10
}
)

spify_parser_params.append(
{
    'name':'thistEcho',
    'type':'bool',
    'longDesc' : "Echo time history to stdout (1=true, 0=false)",
    'defaultValue' : 1
}
)

spify_parser_params.append(
{
    'name':'perturbationTimeHistory',
    'type':'bool',
    'shortDesc' : "Write out time history for all perturbations",
    'defaultValue' : 0
}
)

spify_parser_params.append(
{
    'name':'taskProgressPeriod',
    'type':'int',
    'longDesc' : "The period at which the task progress is reported to stdout",
    'defaultValue' : 50
}
)

spify_parser_params.append(
{
    'name':'outFile',
    'type':'string',
    'longDesc' : "Output file containing the time to max dP/dt for each set of A-Factor multipliers"
}
)

spify_parser_params.append(
{
    'name':'checkFile',
    'type':'string',
    'shortDesc' : "Output file with the A-Factor perturbation stats by reaction"
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
    'defaultValue' : 'MOLE_FRACTION'
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
    'name':"maxKrylovDim",
    'type':'int',
    'shortDesc' : "Maximum size of the Krylov subspace",
    'defaultValue' : 5
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

spify_parser_params.append(
{
    'name':"hrrIntegralPower",
    'type':'v_double',
    'shortDesc' : "Heat release rate integral powers",
    'defaultValue' : [ 2.0, 4.0, 8.0, 16.0 ]
}
)

spify_parser_params.append(
{
    'name':"quadrature_atol",
    'type':'double',
    'shortDesc' : "Absolute tolerance for the heat release quadrature variables",
    'defaultValue' : 1.0e-16
}
)
spify_parser_params.append(
{
    'name':"quadrature_rtol",
    'type':'double',
    'shortDesc' : "Relative tolerance for the heat release quadrature variables",
    'defaultValue' : 1.0e-6
}
)

spify_parser_params.append(
{
    'name':"quadrature_controls_step",
    'type':'bool',
    'shortDesc' : "Use the quadrature variable error to control the timestep",
    'defaultValue' : 1
}
)

#Generate parser code
spg().generate(spify_parser_name,spify_parser_params)


#Done.
