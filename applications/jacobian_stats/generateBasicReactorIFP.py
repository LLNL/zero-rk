#!/usr/bin/env python

import os

from SpifyParserGenerator import SpifyParserGenerator as spg

#Name your parser
spify_parser_name = "BasicReactorIFP"

spify_parser_params = []

#Specify parameters
spify_parser_params.append(
{
    'name':'mechFile',
    'type':'string',
    'shortDesc' : "Chemkin Format Mechansim File [input]"
}
)

spify_parser_params.append(
{
    'name':'thermFile',
    'type':'string',
    'shortDesc' : "Chemkin Format Thermodynamics File [input]"
}
)

spify_parser_params.append(
{
    'name':'mechLogFile',
    'type':'string',
    'shortDesc' : "Mechanism Parser Log File [output]",
    'defaultValue' : os.devnull
}
)

spify_parser_params.append(
{
    'name':'mechStatFile',
    'type':'string',
    'shortDesc' : "Mechanism Statistics File [output]",
}
)

spify_parser_params.append(
{
    'name':'jacobianRawFile',
    'type':'string',
    'shortDesc' : "File containing the raw jacobian data [output]",
}
)

spify_parser_params.append(
{
    'name':'jacobianStatFile',
    'type':'string',
    'shortDesc' : "File containing the jacobian element distributions [output]",
}
)

spify_parser_params.append(
{
    'name':'integratorProbFile',
    'type':'string',
    'shortDesc' : "File containing additional information about CVode errors and warnings [output]",
}
)

spify_parser_params.append(
{
    'name':'outFile',
    'type':'string',
    'shortDesc' : "Output file containing the time history"
}
)

spify_parser_params.append(
{
    'name':'computeEigenvalues',
    'type':'bool',
    'longDesc' : "Compute the eigenvalues of every numerical Jacobian M = I - gamma*J, where J is the ODE Jacobian",
    'defaultValue' : 0  
}
)

spify_parser_params.append(
{
    'name':'eigenvalueFile',
    'type':'string',
    'longDesc' : "File name to store the eigenvalue stats if computeEigenvalues is set to y",
    'defaultValue' : os.devnull
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
    'name':'maxTime',
    'type':'double',
    'shortDesc' : "Maximum integration time"
}
)
spify_parser_params.append(
{
    'name':'printDeltaTime',
    'type':'double',
    'shortDesc' : "Print the current state every printDeltaTime [s]"
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
    'name':"use_unimolecular_limit",
    'type':'bool',
    'longDesc' : "Use [y|n] the value specified as the unimolecular_limit to prevent the unimolecular rate coefficient exceeding a certain value",
    'defaultValue' : 0
}
)
spify_parser_params.append(
{
    'name':"unimolecular_limit",
    'type':'double',
    'longDesc' : "Value L [Hz] that the unimolecular rate coefficient (K_uni) can not exceed  K_limit = K_uni*(L/(K_uni+L))",
    'defaultValue' : 1.0e16
}
)

spify_parser_params.append(
{
    'name':"use_bimolecular_limit",
    'type':'bool',
    'longDesc' : "Use [y|n] the value specified as the bimolecular_limit to prevent the bimolecular rate from exceeding the collision rate by a certain factor",
    'defaultValue' : 0
}
)
spify_parser_params.append(
{
    'name':"bimolecular_limit",
    'type':'double',
    'longDesc' : "Value L [-] that the bimolecular rate coefficient normalized by the bimolecular collision rate  K_bi/K_coll) can not exceed  K_limit = K_bi*(L/(K_bi/K_coll+L))",
    'defaultValue' : 100.0
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

spify_parser_params.append(
{
    'name':"scan_temperature_ref",
    'type':'double',
    'longDesc' : "Reference temperature [K] at which hard sphere diameters are estimated for computing binary step reaction probabilities",
    'defaultValue' : 1000.0  
}
)

spify_parser_params.append(
{
    'name':"min_scan_temperature",
    'type':'double',
    'shortDesc' : "Minimum temperature [K] to start scan for a mechanism check",
    'defaultValue' : 300.0  
}
)

spify_parser_params.append(
{
    'name':"max_scan_temperature",
    'type':'double',
    'shortDesc' : "Maximum temperature [K] to end scan for a mechanism check",
    'defaultValue' : 3000.0  
}
)

spify_parser_params.append(
{
    'name':"num_scans",
    'type':'int',
    'shortDesc' : "Number of temperatures in scan for mechanism check",
    'defaultValue' : 28  
}
)

spify_parser_params.append(
{
    'name':"min_unimolecular_rate",
    'type':'double',
    'shortDesc' : "Minimum unimolecular rate [Hz] to report",
    'defaultValue' : 1.0e14
}
)

spify_parser_params.append(
{
    'name':"min_bimolecular_prob",
    'type':'double',
    'shortDesc' : "Minimum bimolecular probability [-] to report",
    'defaultValue' : 1.0
}
)

spify_parser_params.append(
{
    'name':"min_jacobian_raw",
    'type':'double',
    'longDesc' : "Minimum species Jacobian term [Hz] to report in jacobianRawFile",
    'defaultValue' : 1.0e15
}
)

spify_parser_params.append(
{
    'name':"min_df_dtemp_raw",
    'type':'double',
    'longDesc' : "Minimum d*/dTemp Jacobian term [Hz] to report in jacobianRawFile",
    'defaultValue' : 1.0e15
}
)


spify_parser_params.append(
{
    'name':"min_cvode_species_error",
    'type':'double',
    'longDesc' : "Minimum cvode scaled error of the limiting species above which species statistics are collected.  A value of scaled error of one is considered to be at the relative and absolute error tolerances",
    'defaultValue' : 0.25
}
)

spify_parser_params.append(
{
    'name':"min_cvode_reaction_error",
    'type':'double',
    'longDesc' : "Minimum cvode scaled error of the limiting reaction above which species statistics are collected.  A value of scaled error of one is considered to be at the relative and absolute error tolerances",
    'defaultValue' : 0.25
}
)

spify_parser_params.append(
{
    'name':"pollute_zero_mole_fractions",
    'type':'bool',
    'longDesc' : "Initialize state array with small non-zero mole fractions for all species",
    'defaultValue' : 0
}
)




#Generate parser code (will be created in build directory)
spg().generate(spify_parser_name,spify_parser_params)

# Make master file (will be created in build directory)
spg().make_master_file(spify_parser_name,spify_parser_params)

#Done.
