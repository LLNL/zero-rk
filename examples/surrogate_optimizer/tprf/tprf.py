#!/bin/env python


from zerork.config import ZERORK_ROOT
from zerork.surrogate_optimizer import *


#TODO: These are not distributable yet and so are not in the repo
mech_file=os.path.join(ZERORK_ROOT,'share','zerork','mechanisms','gasoline_surrogates','DofF.inp')
therm_file=os.path.join(ZERORK_ROOT,'share','zerork','mechanisms','gasoline_surrogates','DofF.dat')

db_file='./Database.csv'
mixture = SurrogateMixture(db_file,mech_file,therm_file)

#Define Targets:
targets = []

# RON
RON_target=96
RON_weight=2.0
targets.append(RON_NN_Target(RON_target,RON_weight))

# MON
MON_target=88
MON_weight=2.0
targets.append(MON_NN_Target(MON_target,MON_weight))


# Set initial guess
composition = {
'IC8': 0.8,
'NC7H16': 0.15,
'C6H5CH3': 0.05,
}

mixture.optimize(composition,targets,jobname="TPRF",mpi_procs=4)

