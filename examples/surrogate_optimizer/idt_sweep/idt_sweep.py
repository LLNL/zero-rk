#!/bin/env python

import os
import numpy
from tqdm import tqdm
from zerork import SurrogateMixture
from zerork.config import ZERORK_ROOT
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor

def calc_idts(initial_condition):
    """Calculate idts for comp for temperatures and pressures."""
    label, comp, temperatures, pressures, phis, egrs = initial_condition

    stopping_time = 2.0
    delta_temp_ignition = 400.0  #idt temperature rise metric
    oxidizer_mole_fractions = {"O2":0.1, "AR":0.9}
    trace_mole_fractions = {"NO":150.0e-6}
    for sp in comp.keys():
        if sp not in trace_mole_fractions:
            trace_mole_fractions[sp] = 0.0

    idts = mixture.get_idts(comp, temps=temperatures, press=pressures, jobname=run_name,
                            stdout=False,logout=False, phis=phis, egrs=egrs, 
                            stopping_time=stopping_time,
                            delta_temp_ignition=delta_temp_ignition,
                            oxidizer_mole_fractions=oxidizer_mole_fractions,
                            trace_mole_fractions=trace_mole_fractions)

    calc_species = mixture.calc_species
    species_init_mole_fracs = mixture.species_init_mole_fracs
    return (idts, calc_species, species_init_mole_fracs)

def array_output(output):
    with open('{}.dat'.format(run_name),'w') as outfile:
        # Header
        outfile.write("# {:>15s}".format("label"))
        for spn in output[0][1][0]:
            outfile.write(" {:>15s}".format(spn))
        outfile.write(" {:>15s} {:>15s} {:>15s} {:>15s} {:>15s}\n".format("temperature", "pressure","phi","egr", "idt"))

        # Output
        for ic,op in zip(ic_array,output):
            label = ic[0]
            mix = ic[1]
            temp = ic[2]
            press = ic[3]
            phi = ic[4]
            egr = ic[5]
            idts = numpy.reshape(op[0],(len(temp),len(press),len(phi),len(egr)),'F')
            species_init_mole_fracs = numpy.reshape(op[2],(len(temp),len(press),len(phi),len(egr),len(op[1][0])),'F')
            for i, t in enumerate(temp):
                for j, p in enumerate(press):
                    for k, f in enumerate(phi):
                        for l, e in enumerate(egr):
                            outfile.write("  {:>15s}".format(label))
                            for imf in species_init_mole_fracs[i][j][k][l][:]:
                                outfile.write(" {:>15g}".format(imf))
                            idt = idts[i][j][k][l]
                            outfile.write(" {:>15g} {:>15g} {:>15g} {:>15g} {:>15g}\n".format(t,p,f,e,idt))

def matrix_output(output):
    #Make an empty database
    data_dict = {}
    for label in labels:
        data_dict[label] = {}
        for temperature in temp_array:
            data_dict[label][temperature] = {}
            for pressure in press_array:
                data_dict[label][temperature][pressure] = {}
                for phi in phi_array:
                    data_dict[label][temperature][pressure][phi] = {}
                    for egr in egr_array:
                        data_dict[label][temperature][pressure][egr] = None
 

    #Fill the database
    for ic,op in zip(ic_array,output):
        label = ic[0]
        mix = ic[1]
        temp = ic[2]
        press = ic[3]
        phi = ic[4]
        egr = ic[5]
        idts = numpy.reshape(op[0],(len(temp),len(press),len(phi),len(egr)),'F')
        species_init_mole_fracs = numpy.reshape(op[2],(len(temp),len(press),len(phi),len(egr),len(op[1][0])),'F')
        for i, t in enumerate(temp):
            for j, p in enumerate(press):
                for k, f in enumerate(phi):
                    for l, e in enumerate(egr):
                        data_dict[label][t][p][f][e] = idts[i][j][k][l]

    # Output
    with open('{}_matrix.csv'.format(run_name),'w') as outfile:
        for label,mix in zip(labels,mixes):
            for phi in phi_array:
                for egr in egr_array:
                    outfile.write("\n#{} phi={} egr={}: ".format(label,phi,egr))
                    for spn in species_names:
                        outfile.write(" {}{{{}}}".format(spn, mix.get(spn,0.0)))
                    outfile.write("\n{:15},".format(" "))
                    for pressure in press_array:
                        outfile.write("{:15g},".format(pressure))
                    outfile.write("\n")
                    for temp in temp_array:
                        outfile.write("{:15g},".format(temp))
                        for press in press_array:
                            idt = data_dict[label][temp][press][phi][egr]
                            outfile.write("{:15g},".format(idt))
                        outfile.write("\n")

# Only do postprocessing/output on rank 0
if __name__ == "__main__":
    comm=MPI.COMM_WORLD
    rank=comm.Get_rank()

    # Make mixture
    mech_file=os.path.join(ZERORK_ROOT,'share','zerork','mechanisms','gasoline_surrogates','DofF.inp')
    therm_file=os.path.join(ZERORK_ROOT,'share','zerork','mechanisms','gasoline_surrogates','DofF.dat')
    db_file='../tprf/Database.csv'
    run_name = 'idt_sweep'
    mixture = SurrogateMixture(db_file,mech_file,therm_file)

    if rank == 0:
        #Conditions
        temp_array = numpy.arange(700.0, 1225.0, 100.0)
        press_array = numpy.arange(10.0*101325, 55.0*101325, 10.0*101325)
        phi_array = numpy.arange(0.5,1.5,0.5)
        egr_array = numpy.arange(0.0,0.3,0.1)

        labels = []
        mixes = []
        species_names = []
        # CSV input for mixtures
        csv_filename = './mixtures.csv'
        with open(csv_filename,'r') as csv_file:
            header = csv_file.readline().strip()
            species_names = header.split(',')[1:]
            for line in csv_file:
                values = line.split(',')
                if len(values) != len(species_names) + 1:
                    print("Incorrect number of values in CSV file (header length = {0}, line length = {1}".format(len(species_names), len(values)))
                labels.append(values[0])
                mix = {}
                for sp, vf in zip(species_names, values[1:]):
                    mix[sp] = float(vf)/100.0 #N.B. file values are percents, but we want fractional
                mixes.append(mix)


        ic_array = []
        for l,m in zip(labels,mixes):
            #Sub-arrays are to allow for paralellization and load-balancing
            for t in temp_array:
                for p in press_array:
                    ic_array.append([l,m,[t],[p],phi_array,egr_array])
    else:
        ic_array = None

    output = None
    with MPICommExecutor(MPI.COMM_WORLD, root=0) as executor:
        if executor is not None:
            output = list(tqdm(executor.map(calc_idts,ic_array), total=len(ic_array)))

    if rank == 0:
        array_output(output)
        matrix_output(output)


