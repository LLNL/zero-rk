#!/usr/bin/env python

from SpifyParserGenerator import SpifyParserGenerator as spg

#Name your parser
spify_parser_name = "ZeroRKCFDPluginTesterIFP"

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
    'defaultValue' : '/dev/null'
}
)

spify_parser_params.append(
{
    'name':'outFile',
    'type':'string',
    'shortDesc' : "Ignition Delay Output File",
}
)

spify_parser_params.append(
{
    'name':'nReactors',
    'type':'int',
    'shortDesc' : "Number of reactors to solve",
    'boundMin': 1
}
)

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
    'name':'deltaTign',
    'type':'double',
    'defaultValue': 400.0,
    'boundMin': 50.0,
    'boundMax': 100000.0
}
)

spify_parser_params.append(
{
    'name':'reactorTime',
    'type':'double',
    'shortDesc' : "Total integration time"
}
)

spify_parser_params.append(
{
    'name':'printDt',
    'type':'double',
    'shortDesc' : "Print Interval"
}
)

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
    'name':"refTemp",
    'type':'double',
    'shortDesc' : "Ignition Delat Metric: Maximum species concentration",
    'defaultValue' : 1000.0
}
)

spify_parser_params.append(
{
    'name':"TMin",
    'type':'double',
    'shortDesc' : "Array of initial temperatures to sweep over",
}
)

spify_parser_params.append(
{
    'name':"TMax",
    'type':'double',
    'shortDesc' : "Array of initial temperatures to sweep over"
}
)

spify_parser_params.append(
{
    'name':"pMin",
    'type':'double',
    'shortDesc' : "Array of initial temperatures to sweep over"
}
)

spify_parser_params.append(
{
    'name':"pMax",
    'type':'double',
    'shortDesc' : "Array of initial temperatures to sweep over"
}
)

spify_parser_params.append(
{
    'name':"phiMin",
    'type':'double',
    'shortDesc' : "Array of initial temperatures to sweep over"
}
)

spify_parser_params.append(
{
    'name':"phiMax",
    'type':'double',
    'shortDesc' : "Array of initial temperatures to sweep over"
}
)

spify_parser_params.append(
{
    'name':"egrMin",
    'type':'double',
    'shortDesc' : "Array of initial temperatures to sweep over"
}
)

spify_parser_params.append(
{
    'name':"egrMax",
    'type':'double',
    'shortDesc' : "Array of initial temperatures to sweep over"
}
)

spify_parser_params.append(
{
    'name':"precThresh",
    'type':'double',
    'shortDesc' : "Preconditioner threshold values"
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
    'name':"oneStepMode",
    'type':'bool',
    'shortDesc' : "Switch: integrate in one step mode.",
    'defaultValue' : 0
}
)

spify_parser_params.append(
{
    'name':"printAllSteps",
    'type':'bool',
    'shortDesc' : "Switch: print every integrator step.",
    'defaultValue' : 0
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
    'name':"partialPivotThresh",
    'type':'double',
    'shortDesc' : "partial pivot threshold.",
    'defaultValue' : 0.0
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
    'name':"maxGammaOrderChange",
    'type':'double',
    'shortDesc' : "??",
    'defaultValue' : 3.0
}
)

spify_parser_params.append(
{
    'name':"cvStabLimDet",
    'type':'bool',
    'shortDesc' : "??",
    'defaultValue' : 0
}
)

spify_parser_params.append(
{
    'name':"maxNumNonLinIters",
    'type':'int',
    'shortDesc' : "??",
    'defaultValue' : 3
}
)

spify_parser_params.append(
{
    'name':"maxOrd",
    'type':'int',
    'shortDesc' : "??",
    'defaultValue' : 5
}
)

spify_parser_params.append(
{
    'name':"cvEpsLin",
    'type':'double',
    'shortDesc' : "??",
    'defaultValue' : 0.1
}
)

spify_parser_params.append(
{
    'name':"cvNlConvCoeff",
    'type':'double',
    'shortDesc' : "??",
    'defaultValue' : 0.05
}
)

spify_parser_params.append(
{
    'name':"reInitOnPrint",
    'type':'bool',
    'shortDesc' : "Re-initialize CVODE at every print interval",
    'defaultValue' : 0
}
)

spify_parser_params.append(
{
    'name':"dumpStateOnPrint",
    'type':'bool',
    'shortDesc' : "Dump system state at every print interval",
    'defaultValue' : 0
}
)

spify_parser_params.append(
{
    'name':"stateFiles",
    'type':'v_string',
    'shortDesc' : "Files to read initial states from.",
    'defaultValue' : []   #empty
}
)

spify_parser_params.append(
{
    'name':"linearSolver",
    'type':'string',
    'shortDesc' : "",
    'defaultValue' : "IterativeSparse",
    'discreteValues': ["IterativeSparse","DirectDense","DirectDenseDVD"]
}
)

#Generate parser code
spg().generate(spify_parser_name,spify_parser_params)
#spg().make_master_file(spify_parser_name,spify_parser_params)

#Done.
