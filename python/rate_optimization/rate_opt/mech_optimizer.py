
import os
import shutil
import tempfile
import numpy as np
from multiprocessing import Pool

from collections import OrderedDict

def parallel_opt_fn(args):
    app, mf, tf = args
    return app.opt_fn(mf, tf)

class mech_optimizer(object):
    def __init__(self, mech_file, therm_file, nprocs=1):
        self.verbose = False
        self.printed_header = False
        self.prelines = []
        self.postlines = []
        self.reactions = []
        self.reaction_eqns = []
        self.opt_apps = []
        self.nprocs = nprocs
        self.duplicates = dict()
        self.opt_bounds = []
        self.param_fns = []
        self.opt_rxn_idxs = set()
        self.mech_file=None #we don't run with original mech at all
        self.therm_file = therm_file
        self.parse(mech_file)
    def parse(self, mech_file):
        with open(mech_file,'r', encoding='utf8') as mf:
            self.lines = mf.readlines()
        reaction_section=False
        post_reaction_section=False
        for line in self.lines:
            comment = "!".join(line.split("!")[1:])
            linex = line.split("!")[0].strip()
            if 'REAC' == line.strip().upper()[0:4]:
                reaction_section=True
                self.prelines.append(line)
#                for token in line.strip().upper().split():
#                    if token in self.e_unit_conversion_map:
#                        self.e_unit_conversion = self.e_unit_conversion_map(token)
#                        if self.e_unit_conversion is None:
#                            raise ValueError("Unsupported energy units.")
                continue
            if reaction_section:
                if 'END' == line.strip().upper()[0:3]:
                    reaction_section=False
                    post_reaction_section=True
                    self.postlines.append(line)
                    continue
                if "=" in linex:
                    tokens=linex.split()
                    A, n, Ea = map(float,tokens[-3:])
                    eqn = ''.join(tokens[:-3])
                    rxn = {'eqn': eqn, 'A': A, 'n': n, 'Ea': Ea, 'aux': []}
                    rxn.update({'Aorig':A, 'Eaorig':Ea})
                    rxn['comment'] = comment
                    rxn['aux_comments'] = []
                    rxn['line'] = line
                    if eqn in self.reaction_eqns:
                        idxs = [i for i, x in enumerate(self.reaction_eqns) if x == eqn]
                        self.duplicates[eqn] = idxs + [len(self.reactions)]
                    self.reaction_eqns.append(rxn['eqn'])
                    self.reactions.append(rxn)
                else:
                    self.reactions[-1]['aux'].append(line)
                    self.reactions[-1]['aux_comments'].append(comment)
                    if 'PLOG' == linex.upper()[0:4]:
                        plog_vals = list(map(float,linex.split("/")[1].split()))
                        if 'plog_orig' in self.reactions[-1]:
                            self.reactions[-1]['plog_orig'].append(list(plog_vals))
                            self.reactions[-1]['plog_mod'].append(list(plog_vals))
                        else:
                            self.reactions[-1]['plog_orig'] = [list(plog_vals)]
                            self.reactions[-1]['plog_mod'] = [list(plog_vals)]
            elif post_reaction_section:
                self.postlines.append(line)
            else:
                self.prelines.append(line)
    def mult_A_fn(self, rxn_idx_str):
        def fn(A_mult):
            for ridx in map(int, rxn_idx_str.split("-")):
                self.reactions[ridx]['A'] = self.reactions[ridx]['Aorig']*A_mult
                if 'plog_mod' in self.reactions[ridx]:
                    for lidx in range(len(self.reactions[ridx]['plog_mod'])):
                        self.reactions[ridx]['plog_mod'][lidx][1] = \
                              self.reactions[ridx]['plog_orig'][lidx][1]*A_mult
        return fn
    def add_Ea_fn(self, rxn_idx_str):
        def fn(Ea_delta):
            for ridx in map(int, rxn_idx_str.split("-")):
                self.reactions[ridx]['Ea'] = self.reactions[ridx]['Eaorig']+Ea_delta
                if 'plog_mod' in self.reactions[ridx]:
                    for lidx in range(len(self.reactions[ridx]['plog_mod'])):
                        self.reactions[ridx]['plog_mod'][lidx][3] = \
                              self.reactions[ridx]['plog_orig'][lidx][3]+Ea_delta
        return fn
    def set_rxn_opt(self, rxn_eqn, A_mult=2, Ea_delta=None):
        rxn_eqn = ''.join(rxn_eqn.strip().split()) #Remove whitespace
        if(rxn_eqn not in self.reaction_eqns):
            raise ValueError(f"Unrecognized reaction: {rxn_eqn}")
        elif(A_mult is None and Ea_delta is None):
            raise ValueError(f"Please set a value for either A_mult or Ea_delta: {rxn_eqn}")
        else:
            if(rxn_eqn in self.duplicates):
                idxs = self.duplicates[rxn_eqn]
                rxn_idx_str = "-".join(map(str,idxs))
            else:
                rxn_idx_str = str(self.reaction_eqns.index(rxn_eqn))
            for ridx in map(int, rxn_idx_str.split("-")):
                if ridx in self.opt_rxn_idxs:
                    print(f"Reaction already included ({rxn_eqn}). Ignoring.")
                    return
                self.opt_rxn_idxs.add(ridx)
            if A_mult is not None:
                self.opt_bounds.append((1/A_mult, A_mult))
                self.param_fns.append( self.mult_A_fn(rxn_idx_str) )
            if Ea_delta is not None:
                self.opt_bounds.append((-Ea_delta, Ea_delta))
                self.param_fns.append( self.add_Ea_fn(rxn_idx_str) )
    def get_opt_bounds(self):
        return self.opt_bounds
    def get_num_opt_vars(self):
        return len(self.opt_bounds)
    def add_opt_app(self, app):
        self.opt_apps.append(app)
    def set_verbose(self, v):
        self.verbose = bool(v)
    def write_mech(self, params, mech_file):
        for p, fn in zip(params, self.param_fns):
            fn(p)
        with open(mech_file,'w') as mfile:
           for line in self.prelines:
               mfile.write(line)
           for idx, rxn in enumerate(self.reactions):
               if idx in self.opt_rxn_idxs:
                   mfile.write(f"{rxn['eqn']} {rxn['A']:0.4e} {rxn['n']:0.6g} {rxn['Ea']:0.6g}")
                   mfile.write(f"   ! Zero-RK optimized from: {rxn['Aorig']:0.4e} {rxn['n']:0.6g} {rxn['Eaorig']:0.6g}")
                   if len(rxn['comment']) > 0:
                       mfile.write(" !" + rxn['comment'])
                   else:
                       mfile.write("\n")
               else:
                   mfile.write(rxn['line'])
               lidx = 0
               for l,c in zip(rxn['aux'],rxn['aux_comments']):
                   if 'PLOG' == l.strip().upper()[0:4] and idx in self.opt_rxn_idxs:
                       p, A, n, Ea = self.reactions[idx]['plog_mod'][lidx]
                       ll = f"PLOG / {p:4.2e} {A:0.4e} {n:0.6g} {Ea:0.6g} /"
                       p, A, n, Ea = self.reactions[idx]['plog_orig'][lidx]
                       ll += f" ! Zero-RK optimized from {A:0.4e} {n:0.6g} {Ea:0.6g}"
                       if(len(c)>0):
                           ll += " !" + c
                       else:
                           ll += "\n"
                       mfile.write(ll)
                       lidx += 1
                   else: 
                       mfile.write(l)
           for line in self.postlines:
               mfile.write(line)
    def opt_fn(self, params):
        tmpdir = tempfile.mkdtemp(dir='.')
        mech_file=os.path.join(tmpdir,'chem.inp')
        self.write_mech(params,mech_file)
        val = 0
        try:
            if self.nprocs > 1:
                inputs = [(app, mech_file, self.therm_file) for app in self.opt_apps]
                with Pool(self.nprocs) as p:
                    values = p.map(parallel_opt_fn, inputs)
                val = np.sum(values)
            else:
                for app in self.opt_apps:
                    val += app.opt_fn(mech_file, self.therm_file)
            val = val / len(self.opt_apps)
        finally:
            shutil.rmtree(tmpdir)
        if(self.verbose):
            if not self.printed_header:
                print("#"+",".join([f"p{i}" for i,p in enumerate(params)])+ f",objective")
                self.printed_header = True
            print(",".join([f"{p:5.6f}" for p in params])+ f",{val:0.8g}")
        return val



