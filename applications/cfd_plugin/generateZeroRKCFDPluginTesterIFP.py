#!/usr/bin/env python

import os,sys

spify_parser_name = "ZeroRKCFDPluginTesterIFP"

spify_parser_params = []
spify_parser_params.append(
{
    'name':'mechanism_file',
    'type':'string',
    'shortDesc' : "Chemkin Format Mechansim File"
}
)

spify_parser_params.append(
{
    'name':'thermo_file',
    'type':'string',
    'shortDesc' : "Chemkin Format Thermodynamics File"
}
)

spify_parser_params.append(
{
    'name':'zerork_cfd_plugin_input',
    'type':'string',
    'shortDesc' : "Zero-RK Configuration File"
}
)


spify_parser_params.append(
{
    'name':'n_reactors',
    'type':'int',
    'shortDesc' : "Number of reactors to solve",
    'boundMin': 1
}
)

spify_parser_params.append(
{
    'name':'n_print_reactors',
    'type':'int',
    'shortDesc' : "Number of reactors to print history",
    'defaultValue': 0,
    'boundMin': 0
}
)

spify_parser_params.append(
{
    'name':'fuel_composition',
    'type':'m_string_double',
    'shortDesc' : "Fuel composition map"
}
)

spify_parser_params.append(
{
    'name':'oxidizer_composition',
    'type':'m_string_double',
    'shortDesc' : "Oxidizer composition map"
}
)

spify_parser_params.append(
{
    'name':'solution_time',
    'type':'double',
    'shortDesc' : "Total integration time"
}
)

spify_parser_params.append(
{
    'name':'print_dt',
    'type':'double',
    'shortDesc' : "Print Interval"
}
)

spify_parser_params.append(
{
    'name':'n_steps',
    'type':'int',
    'shortDesc' : "Number of 'cfd' steps",
    'defaultValue': 1,
    'boundMin': 1
}
)

spify_parser_params.append(
{
    'name':'delta_temperature_ignition',
    'type':'double',
    'defaultValue': 0.0,
    'boundMin': 0.0,
    'boundMax': 100000.0
}
)

spify_parser_params.append(
{
    'name':"temperature_min",
    'type':'double',
    'shortDesc' : "Array of initial temperatures to sweep over",
}
)

spify_parser_params.append(
{
    'name':"temperature_max",
    'type':'double',
    'shortDesc' : "Array of initial temperatures to sweep over"
}
)

spify_parser_params.append(
{
    'name':"pressure_min",
    'type':'double',
    'shortDesc' : "Array of initial temperatures to sweep over"
}
)

spify_parser_params.append(
{
    'name':"pressure_max",
    'type':'double',
    'shortDesc' : "Array of initial temperatures to sweep over"
}
)

spify_parser_params.append(
{
    'name':"phi_min",
    'type':'double',
    'shortDesc' : "Array of initial temperatures to sweep over"
}
)

spify_parser_params.append(
{
    'name':"phi_max",
    'type':'double',
    'shortDesc' : "Array of initial temperatures to sweep over"
}
)

spify_parser_params.append(
{
    'name':"egr_min",
    'type':'double',
    'shortDesc' : "Array of initial temperatures to sweep over"
}
)

spify_parser_params.append(
{
    'name':"egr_max",
    'type':'double',
    'shortDesc' : "Array of initial temperatures to sweep over"
}
)

spify_parser_params.append(
{
    'name':"state_files",
    'type':'v_string',
    'shortDesc' : "Files to read initial states from.",
    'defaultValue' : []   #empty
}
)

spify_parser_params.append(
{
    'name':"state_files_cfd",
    'type':'v_string',
    'shortDesc' : "CFD files to read initial states from.",
    'defaultValue' : []   #empty
}
)

spify_parser_params.append(
{
    'name':"reactor_history_file_prefix",
    'type':'string',
    'shortDesc' : "Filename prefix for history files.",
    'defaultValue' : "reactor_history"
}
)

spify_parser_params.append(
{
    'name':"log_species",
    'type':'v_string',
    'shortDesc' : "Species names to includie in history (defaults to fuel and oxidizer species.",
    'defaultValue' : []   #empty
}
)

spify_parser_params.append(
{
    'name':"app_owns_aux_fields",
    'type':'bool',
    'shortDesc' : "If true, application over-rides library storage for auxilliary fields (reactor cost, gpu solve, and dP/dt).",
    'defaultValue' : 0
}
)

spify_parser_params.append(
{
    'name':"dpdt",
    'type':'double',
    'shortDesc' : "Pressure gradient in Pa/s",
    'defaultValue' : 0.0
}
)

spify_parser_params.append(
{
    'name':"e_src",
    'type':'double',
    'shortDesc' : "Energy source in J/kg/s",
    'defaultValue' : 0.0
}
)

spify_parser_params.append(
{
    'name':"y_src",
    'type':'double',
    'shortDesc' : "Species mass fraction source in 1/s",
    'defaultValue' : 0.0
}
)

spify_parser_params.append(
{
    'name':"constant_volume",
    'type':'int',
    'shortDesc' : "Energy form: constant volume or constant pressure.",
    'defaultValue' : 1
}
)

spify_parser_params.append(
{
    'name':"stop_after_ignition",
    'type':'int',
    'shortDesc' : "Stop integration once ignition criteria is reached.",
    'defaultValue' : 0
}
)

spify_parser_params.append(
{
    'name':"batched",
    'type':'int',
    'shortDesc' : "Solve multiple reactors in one call to plugin.",
    'defaultValue' : 1
}
)

from SpifyParserGenerator import SpifyParserGenerator as spg

spg().generate(spify_parser_name,spify_parser_params)
#spg().make_master_file(spify_parser_name,spify_parser_params)


