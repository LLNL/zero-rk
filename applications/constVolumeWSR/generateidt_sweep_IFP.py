#!/usr/bin/env python

import os,sys

#Root of spify src directory
SPIFY_SRC_DIR = '${HOME}/src/spify'

#Name your parser
spify_parser_name = "idt_sweep_IFP"

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
    'name':'logFile',
    'type':'string',
    'shortDesc' : "Mechanism Parser Log File",
    'defaultValue' : '/dev/null'
}
)

spify_parser_params.append(
{
    'name':'idtFile',
    'type':'string',
    'shortDesc' : "Ignition Delay Output File",
}
)

spify_parser_params.append(
{
    'name':'thistFile',
    'type':'string',
    'shortDesc' : "Temperature History Output File",
}
)

spify_parser_params.append(
{
    'name':'fuel_mole_fracs',
    'type':'m_string_double',
    'shortDesc' : "Fuel mole fractions"
}
)

spify_parser_params.append(
{
    'name':'oxidizer_mole_fracs',
    'type':'m_string_double',
    'shortDesc' : "Oxidizer mole fractions"
}
)

spify_parser_params.append(
{
    'name':'trace_mole_fracs',
    'type':'m_string_double',
    'shortDesc' : "Trace mole fractions",
    'defaultValue' : {}
}
)

spify_parser_params.append(
{
    'name':'delta_T_ignition',
    'type':'double',
    'defaultValue': 400.0,
    'boundMin': 50.0,
    'boundMax': 100000.0
}
)

spify_parser_params.append(
{
    'name':'stop_time',
    'type':'double',
    'shortDesc' : "Total integration time"
}
)

spify_parser_params.append(
{
    'name':'print_time',
    'type':'double',
    'shortDesc' : "Print Interval"
}
)

spify_parser_params.append(
{
    'name':'max_internal_dt',
    'type':'double',
    'shortDesc' : "Maximum Internal Integrator Step Size",
    'defaultValue' : 0.05
}
)

spify_parser_params.append(
{
    'name' :'max_internal_steps',
    'type':'int',
    'shortDesc' : "Maximum Number of Integrator Steps",
    'defaultValue':1000000
}
)

spify_parser_params.append(
{
    'name':'relative_tolerance',
    'type':'double',
    'shortDesc':"Relative Integrator Tolerance",
    'defaultValue':1.0e-8
}
)

spify_parser_params.append(
{
    'name':"absolute_tolerance",
    'type':'double',
    'shortDesc':"Absolute Integrator Tolerance",
    'defaultValue':1.0e-20
}
)


spify_parser_params.append(
{
    'name':"reference_temperature",
    'type':'double',
    'shortDesc' : "Ignition Delat Metric: Maximum species concentration",
    'defaultValue' : 1000.0
}
)


spify_parser_params.append(
{
    'name':"max_cvode_fails1",
    'type':'int',
    'shortDesc' : "Max allowed CVODE failures at original preconditioner threshold",
    'defaultValue' : 5
}
)

spify_parser_params.append(
{
    'name':"max_cvode_fails2",
    'type':'int',
    'shortDesc' : "Max allowed CVODE failures at safety preconditioner threshold",
    'defaultValue' : 2
}
)

spify_parser_params.append(
{
    'name':"safety_threshold",
    'type':'double',
    'shortDesc' : "Safety preconditionary threshold (used after max_cvode_fails1 failures)",
    'defaultValue' : 1.0e-6
}
)

spify_parser_params.append(
{
    'name':"continue_after_ignition",
    'type':'int',
    'shortDesc' : "Continue integration to stop time after ignition delay criteria met",
    'defaultValue' : 1
}
)


spify_parser_params.append(
{
    'name':"initial_temperatures",
    'type':'v_double',
    'shortDesc' : "Vector of initial temperatures to use in sweep.",
}
)

spify_parser_params.append(
{
    'name':"initial_pressures",
    'type':'v_double',
    'shortDesc' : "Vector of initial pressures to use in sweep.",
}
)

spify_parser_params.append(
{
    'name':"initial_phis",
    'type':'v_double',
    'shortDesc' : "Vector of initial equivalence ratios to use in sweep.",
}
)

spify_parser_params.append(
{
    'name':"initial_egrs",
    'type':'v_double',
    'shortDesc' : "Vector of initial EGR ratios to use in sweep.",
}
)

spify_parser_params.append(
{
    'name':"preconditioner_thresholds",
    'type':'v_double',
    'shortDesc' : "Vector of preconditioner thresholds to use in sweep.",
    'defaultValue' : [1.0e-3]
}
)

spify_parser_params.append(
{
    'name':"max_krylov_dimensions",
    'type':'v_int',
    'shortDesc' : "Vector of Krylov dimensions to use in sweep.",
    'defaultValue' : [5]
}
)



spify_parser_params.append(
{
    'name':"one_step_mode",
    'type':'int',
    'shortDesc' : "Switch: integrate in one step mode.",
    'defaultValue' : 0
}
)

spify_parser_params.append(
{
    'name':"long_output",
    'type':'int',
    'shortDesc' : "Switch: longer floating point output data.",
    'defaultValue' : 0
}
)

spify_parser_params.append(
{
    'name':"fake_update",
    'type':'int',
    'shortDesc' : "Switch: fakeupdate.",
    'defaultValue' : 0
}
)

spify_parser_params.append(
{
    'name':"incomplete_LU",
    'type':'int',
    'shortDesc' : "Switch: ILU.",
    'defaultValue' : 0
}
)

spify_parser_params.append(
{
    'name':"threshold_type",
    'type':'int',
    'shortDesc' : "Switch: Preconditioner Threshold Type.",
    'defaultValue' : 1
}
)

spify_parser_params.append(
{
    'name':"partial_pivot_threshold",
    'type':'double',
    'shortDesc' : "partial pivot threshold.",
    'defaultValue' : 0.0
}
)

spify_parser_params.append(
{
    'name':"permutation_type",
    'type':'int',
    'shortDesc' : "perm type.",
    'defaultValue' : 1
}
)

spify_parser_params.append(
{
    'name':"eps_lin",
    'type':'double',
    'shortDesc' : "??",
    'defaultValue' : 0.1
}
)

spify_parser_params.append(
{
    'name':"nonlinear_convergence_coefficient",
    'type':'double',
    'shortDesc' : "??",
    'defaultValue' : 0.05
}
)

spify_parser_params.append(
{
    'name':"dump_jacobian",
    'type':'int',
    'defaultValue' : 0
}
)


#Make sure we can import SpifyParserGenerator
sys.path.append(os.path.join(os.path.expandvars(SPIFY_SRC_DIR),'src'))

#Import
from SpifyParserGenerator import SpifyParserGenerator as spg

#Generate parser code
spg().generate(spify_parser_name,spify_parser_params)
spg().make_master_file(spify_parser_name,spify_parser_params)


#Done.
