#!/usr/bin/env python

import os,sys

#Root of spify src directory
SPIFY_SRC_DIR = '../../opt/include/spify'

#Name your parser
spify_parser_name = "zerork_reactorIFP"

spify_parser_params = []

#Specify parameters

spify_parser_params.append(
{
    'name':'cuda',
    'type':'bool',
    'shortDesc' : "Enable CUDA?",
    'defaultValue' : 1
}
)

spify_parser_params.append(
{
    'name':'cudaDevice',
    'type':'int',
    'shortDesc' : "CUDA Device to use for Simulation",
    'defaultValue' : -1 # CUDA library picks for us
}
)


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
    'name':'nReactorsMax',
    'type':'int',
    'shortDesc' : "Maximum number of reactors to solve",
    'boundMin': 1,
    'boundMax': 8192,
    'defaultValue': 8192
}
)

spify_parser_params.append(
{
    'name':'nReactorsMin',
    'type':'int',
    'shortDesc' : "Minimum number of reactors to solve",
    'boundMin': 1,
    'defaultValue': 256
}
)

spify_parser_params.append(
{
    'name':'nMatrixReactors',
    'type':'int',
    'shortDesc' : "Number of reactors to solve",
    'boundMin': 1,
    'defaultValue': 8192000 #Bigger than any reasonable number of total reactors
}
)

spify_parser_params.append(
{
    'name':'nMatrixThreads',
    'type':'int',
    'shortDesc' : "Number of threads to use in matrix factor",
    'boundMin': 1,
    'boundMax': 32,
    'defaultValue': 1
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
    'name':'reactorType',
    'type':'string',
    'defaultValue': "CV",
    'discreteValues': ["CV","CP"]
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
    'name':"constantPressure",
    'type':'bool',
    'shortDesc' : "Constant Pressure instead of constant volume",
    'defaultValue' : 0
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
    'shortDesc' : "CVODE Nonlinear solver",
    'defaultValue' : "IterativeSparse",
    'discreteValues': ["IterativeSparse","DirectDense","DirectDenseDVD","DirectSparse","DirectSparseDVD",
                       "IterativeSparseCusolverRf","IterativeSparseCusolverRfBlocked","IterativeSparseCusolverRfBatched",
                       "IterativeSparseCusolverSpBatched"]
}
)

spify_parser_params.append(
{
    'name':"cusolverRf",
    'type':'bool',
    'shortDesc' : "Enable cusolverRf Sparse Solver",
    'defaultValue' : 1
}
)

spify_parser_params.append(
{
    'name':"CUSOLVER_FACTOR_ALG",
    'type':'int',
    'shortDesc' : "cusolverRf Factorization algorithm",
    'defaultValue' : 2,
    'discreteValues': [0,1,2,3]
}
)

spify_parser_params.append(
{
    'name':"CUSOLVER_SOLVE_ALG",
    'type':'int',
    'shortDesc' : "cusolverRf Back Substition algorithm",
    'defaultValue' : 0,
    'discreteValues': [0,1,2]
}
)

spify_parser_params.append(
{
    'name':"cusolverSp",
    'type':'bool',
    'shortDesc' : "Enable cusolverSp Sparse Solver",
    'defaultValue' : 1
}
)

#double max_dt_diff_log = 2.20/(1+3.2*exp(-log_dt_first-2.2))+0.015;
spify_parser_params.append(
{
    'name':"logiA",
    'type':'double',
    'shortDesc' : "logistical function A factor",
    'defaultValue' : 0.015
}
)

spify_parser_params.append(
{
    'name':"logiK",
    'type':'double',
    'shortDesc' : "logistical function K factor",
    'defaultValue' : 2.35
}
)

spify_parser_params.append(
{
    'name':"logiQ",
    'type':'double',
    'shortDesc' : "logistical function Q factor",
    'defaultValue' : 3.2
}
)

spify_parser_params.append(
{
    'name':"logiM",
    'type':'double',
    'shortDesc' : "logistical function M factor",
    'defaultValue' : -2.2
}
)

spify_parser_params.append(
{
    'name':"logiB",
    'type':'double',
    'shortDesc' : "logistical function B factor",
    'defaultValue' : 1.0
}
)

spify_parser_params.append(
{
    'name':"logirNu",
    'type':'double',
    'shortDesc' : "logistical function rNu factor",
    'defaultValue' : 1.0
}
)

spify_parser_params.append(
{
    'name':"min_scaled_dt",
    'type':'double',
    'shortDesc' : "minimum scaled time step for grouping reactors",
    'defaultValue' : 1.0e-5
}
)

#Make sure we can import SpifyParserGenerator
sys.path.append(SPIFY_SRC_DIR)

#Import
from SpifyParserGenerator import SpifyParserGenerator as spg

#Generate parser code
spg().generate(spify_parser_name,spify_parser_params)
#spg().make_master_file(spify_parser_name,spify_parser_params)


#Done.
