#!/usr/bin/env python

import os

from SpifyParserGenerator import SpifyParserGenerator as spg

#Name your parser
spify_parser_name = "AFactorIFP"

spify_parser_params = []

#Specify parameters

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
    'shortDesc' : "Sensitivity output file"
}
)

spify_parser_params.append(
{
    'name':'AFactorMultiplier',
    'type':'double',
    'shortDesc' : "AFactor perturbation multiplier",
    'defaultValue' : 2.0
}
)

spify_parser_params.append(
{
    'name':'doBothDir',
    'type':'bool',
    'shortDesc' : "Switch: compute the AFactor perturbation in both directions",
    'defaultValue' : 0
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
    'name':"stopAfterLastIdtTemp",
    'type':'bool',
    'longDesc' : "Stop the calculation after the last ignition delay time based on the largest change of temperature is found. Note that if this is set to true, then the ignition delay times based on the maximum species mass fractions and the maximum heat release rate may not be correct.",  
    'defaultValue' : 0
}
)

spify_parser_params.append(
{
    'name': 'trackSpeciesMax',
    'type': 'v_string',
    'longDesc': 'Vector of species names to record the time and value at their maximum. If a species is not found, it is ignored.',
    'defaultValue' : [ ]
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
    'type':'v_double',
    'shortDesc' : "Initial temperature [K]",
    'defaultValue' : [1000.0]
}
)

spify_parser_params.append(
{
    'name':"initPres",
    'type':'v_double',
    'shortDesc' : "Initial Pressure [Pa]",
    'defaultValue' : [1.0e5]
}
)

spify_parser_params.append(
{
    'name':"initPhi",
    'type':'v_double',
    'shortDesc' : "Initial equivalence ratio (F/A)/(F/A)_{st} [-]",
    'defaultValue' : [1.0]
}
)

spify_parser_params.append(
{
    'name':"doTPPhiList",
    'type':'bool',
    'longDesc' : "If true, perform the A-Factor sensitivity analysis treating the vectors of initial temperature, pressure and phi as a list. That is, run 0 will be at {T[0], P[0], Phi[0]}; run 1 will be at {T[1], P[1], Phi[1]}; and so forth. If the vectors are of unequal length, then the run list will be truncated to the shortest vector. If false, perform the A-Factor sensitivity as a three dimensional sweep so that the number of runs = (size of the initial tempeture vector) * (size of the initial pressure vector) * (size of the initial phi vector)",
    'defaultValue' : 1
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
