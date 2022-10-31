#!/usr/bin/env python

import os

from SpifyParserGenerator import SpifyParserGenerator as spg

#Name your parser
spify_parser_name = "GSA_AFactorIFP"

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
    'name':'outFile',
    'type':'string',
    'shortDesc' : "Igntion delay output file"
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
    'name':'gsaMatrixFile',
    'type':'string',
    'shortDesc' : "Input file containing the multiplicative A-Factor perturbations"
}
)



"""
spify_parser_params.append(
{
    'name':'nReactors',
    'type':'int',
    'shortDesc' : "Number of reactors to solve",
    'boundMax': 1024
}
)
"""

spify_parser_params.append(
{
    'name':'fuelComp',
    'type':'m_string_double',
    'shortDesc' : "Fuel composition map"
}
)

spify_parser_params.append(
{
    'name':'oxidizerComp',
    'type':'m_string_double',
    'shortDesc' : "Oxidizer composition map"
}
)

spify_parser_params.append(
{
    'name':'reactorType',
    'type':'string',
    'defaultValue': "CV",
    'discreteValues': ["CV","CP"]
}
)

spify_parser_params.append(
{
    'name':'maxTime',
    'type':'double',
    'shortDesc' : "Maximum integration time"
}
)

"""
spify_parser_params.append(
{
    'name':'printDt',
    'type':'double',
    'shortDesc' : "Print Interval"
}
)
"""

spify_parser_params.append(
{
    'name':'maxDtInternal',
    'type':'double',
    'shortDesc' : "Maximum Internal Integrator Step Size",
    'defaultValue' : 0.05
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
    'name':"idtTemps",
    'type':'v_double',
    'shortDesc' : "Ignition Delay Metric: Rise in temperature",  
    'defaultValue' : [ 400.0 ]  #single metric 400K rise
}
)

spify_parser_params.append(
{
    'name':"refTemp",
    'type':'double',
    'shortDesc' : "Reference temperature to normalize the ODE dT/dt equation",
    'defaultValue' : 1000.0  
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
    'name':"initPhi",
    'type':'double',
    'shortDesc' : "Initial equivalence ratio (F/A)/(F/A)_{st} [-]",
}
)
spify_parser_params.append(
{
    'name':"initEgr",
    'type':'double',
    'shortDesc' : "Idealized exhaust gas recirculation fraction [-]",
    'defaultValue' : 0.0   
}
)

spify_parser_params.append(
{
    'name':"precThresh",
    'type':'double',
    'shortDesc' : "Preconditioner threshold value",
    'defaultValue' : 0.0
}
)

spify_parser_params.append(
{
    'name':"krylovDim",
    'type':'int',
    'shortDesc' : "Integrator max krylov dimension",
    'defaultValue' : 5   #CVode default for CVSpgmr
}
)

spify_parser_params.append(
{
    'name':"doFakeUpdate",
    'type':'bool',
    'shortDesc' : "Switch: fakeupdate.",
    'defaultValue' : 0  
}
)

spify_parser_params.append(
{
    'name':"doILU",
    'type':'bool',
    'shortDesc' : "Switch: ILU.",
    'defaultValue' : 0  
}
)

spify_parser_params.append(
{
    'name':"precThreshType",
    'type':'int',
    'shortDesc' : "Switch: Preconditioner Threshold Type.",
    'defaultValue' : 1  
}
)

spify_parser_params.append(
{
    'name':"partialPivotThresh",
    'type':'double',
    'shortDesc' : "partial pivot threshold.",
    'defaultValue' : 0.0  
}
)

spify_parser_params.append(
{
    'name':"permutationType",
    'type':'int',
    'shortDesc' : "perm type.",
    'defaultValue' : 1  
}
)

spify_parser_params.append(
{
    'name':"permThresh",
    'type':'double',
    'shortDesc' : "perm threshold.",
    'defaultValue' : 0.3  
}
)

spify_parser_params.append(
{
    'name':"strictSamePattern",
    'type':'bool',
    'shortDesc' : "Switch",
    'defaultValue' : 0  
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

#Generate parser code
spg().generate(spify_parser_name,spify_parser_params)


#Done.
