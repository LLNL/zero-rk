#!/bin/env python

import numpy
from zerork import *

def get_mixture(db_file="./Database.csv"):
    #We don't need mech or therm unless using IDT for RON/MON
    mech_file=None
    therm_file=None
    mixture = SurrogateMixture(db_file,mech_file,therm_file)
    return mixture

def get_targets():
    #Define Targets:
    targets = []

    CN_target=44
    CN_weight=4.0
    targets.append(CetaneNumberTarget(CN_target,CN_weight))


    Density_target=0.848
    Density_weight=3.0
    targets.append(DensityTarget(Density_target,Density_weight))


    # Distillation curve data
    dc_data_file = "dc_CFA.csv"
    dc_weight = 12.0
    targets.append(DistillationCurveTarget(dc_data_file, dc_weight,units="C", offset=4.9))


    # Carbon Type
    C_data_file = "C_Target.csv"
    C_weight = 0.0
    targets.append(CarbonTypeTarget(C_data_file, C_weight))


    # Hydrogen-Carbon Ratio
    HC_target=1.7825
    HC_weight =3.0
    targets.append(HCTarget(HC_target,HC_weight))

    CT_base_weight = 4.0

    C_data_file = "C_Target.csv"
    C_type_num = 1
    C_weight = CT_base_weight
    targets.append(SingleCarbonTypeTarget(C_data_file, C_type_num, C_weight))

    C_data_file = "C_Target.csv"
    C_type_num = 2
    C_weight = CT_base_weight
    targets.append(SingleCarbonTypeTarget(C_data_file, C_type_num, C_weight))

    C_data_file = "C_Target.csv"
    C_type_num = 3
    C_weight = CT_base_weight/20.0
    targets.append(SingleCarbonTypeTarget(C_data_file, C_type_num, C_weight))

    C_data_file = "C_Target.csv"
    C_type_num = 4
    C_weight = CT_base_weight
    targets.append(SingleCarbonTypeTarget(C_data_file, C_type_num, C_weight))

    C_data_file = "C_Target.csv"
    C_type_num = 5
    C_weight = CT_base_weight/2.0
    targets.append(SingleCarbonTypeTarget(C_data_file, C_type_num, C_weight))

    C_data_file = "C_Target.csv"
    C_type_num = 6
    C_weight = CT_base_weight/2.0
    targets.append(SingleCarbonTypeTarget(C_data_file, C_type_num, C_weight))

    C_data_file = "C_Target.csv"
    C_type_num = 7
    C_weight = CT_base_weight
    targets.append(SingleCarbonTypeTarget(C_data_file, C_type_num, C_weight))

    C_data_file = "C_Target.csv"
    C_type_num = 8
    C_weight = CT_base_weight/4.0
    targets.append(SingleCarbonTypeTarget(C_data_file, C_type_num, C_weight))

    C_data_file = "C_Target.csv"
    C_type_num = 9
    C_weight = CT_base_weight/4.0
    targets.append(SingleCarbonTypeTarget(C_data_file, C_type_num, C_weight))

    C_data_file = "C_Target.csv"
    C_type_num = 10
    C_weight = CT_base_weight/4.0
    targets.append(SingleCarbonTypeTarget(C_data_file, C_type_num, C_weight))

    C_data_file = "C_Target.csv"
    C_type_num = 11
    C_weight = CT_base_weight/8.0
    targets.append(SingleCarbonTypeTarget(C_data_file, C_type_num, C_weight))

    # Limit volume fractions
    vf_max=0.05
    vf_max_weight=10.0
    vf_sp = 'NC18H38'
    targets.append(LimitedVolumeFractionTarget(vf_max,vf_max_weight,vf_sp))

    vf_max=0.05
    vf_max_weight=10.0
    vf_sp = 'NC20H42'
    targets.append(LimitedVolumeFractionTarget(vf_max,vf_max_weight,vf_sp))

    vf_max=0.15
    vf_max_weight=10.0
    vf_sp = 'TDECALIN'
    targets.append(LimitedMassFractionTarget(vf_max,vf_max_weight,vf_sp))

    vf_max=0.15
    vf_max_weight=10.0
    vf_sp = 'A2CH3'
    targets.append(LimitedMassFractionTarget(vf_max,vf_max_weight,vf_sp))

    vf_max=0.10
    vf_max_weight=10.0
    vf_sp = 'TETRA'
    targets.append(LimitedMassFractionTarget(vf_max,vf_max_weight,vf_sp))

    vf_max=0.15
    vf_max_weight=10.0
    vf_sp = 'C6H5C4H9'
    targets.append(LimitedMassFractionTarget(vf_max,vf_max_weight,vf_sp))

    vf_max=0.15
    vf_max_weight=10.0
    vf_sp = 'T124MBZ'
    targets.append(LimitedMassFractionTarget(vf_max,vf_max_weight,vf_sp))

    vf_max=0.15
    vf_max_weight=10.0
    vf_sp = 'NBCH'
    targets.append(LimitedMassFractionTarget(vf_max,vf_max_weight,vf_sp))

    vf_max=0.30
    vf_max_weight=10.0
    vf_sp = 'HMN'
    targets.append(LimitedMassFractionTarget(vf_max,vf_max_weight,vf_sp))

    YSI_target = 122.2
    YSI_weight = 0.1
    targets.append(YSITarget(YSI_target,YSI_weight,mode="CT"))

    return targets


def get_species_list():
    # Set initial guess
    species = [
    'T124MBZ',
    'A2CH3',
    'NC16H34',
    'NC14H30',
    'NC12H26',
    'TETRA',
    'TDECALIN',
    'NC20H42',
    'NC18H38',
    'PBZ'    ,
    'C6H5C4H9',
    'HMN'    ,
    'NBCH'  ,
    ]
    #Use all species in databas
    #species = [sp.name for sp in mixture.species]

    return species

def main():
    run_name = "DieselCFA_ea"
    mixture = get_mixture()
    targets = get_targets()
    species = get_species_list()

    mixture.setTransform("mult")
    init_comp = { sp: 1/len(species) for sp in species }
    mixture.optimize_ea(init_comp,targets,
                        stdout=False,
                        logout=True,
                        jobname="{}".format(run_name))

    print("Evaluating best composition")
    print(f"Distillation prediction will be in {run_name}_best_dist_curve.txt")
    best_comp = { sp.name: vf for sp, vf in zip(mixture.species, mixture.vfs_best[-1][0]) if vf > 0}
    mixture.evaluate(best_comp, targets, jobname=f"{run_name}_best")


if __name__ == "__main__":
    main()

