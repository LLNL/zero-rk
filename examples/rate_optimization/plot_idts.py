#!/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.constrained_layout.use'] = True

#User must create idt_*.dat from sweep.yml

labels = ['orig', 'mod', 'opti']
pretty_labels={"orig": "Original", "mod": "Modified", "opti": "Optimized"}
idts = []
for l in labels:
    try:
        idts.append(np.genfromtxt(f"idt_{l}.dat", comments="#", skip_footer=1))
    except Exception as e:
        print(f"\n\nCouldn't load file idt_{l}.dat. Please generate with sweep.yml\n\n")
        raise e

colors = ['black', 'green', 'blue']
line_styles = ['-', ':', '--']
fig, axes = plt.subplots(1)
for data,l,ls in zip(idts, labels, line_styles):
    pl = pretty_labels[l]
    lw = 3.5
    if l == 'orig':
        lw = 1.5
    phi_array = np.unique(data[:,3])
    P_array = np.unique(data[:,2])
    T_array = np.unique(data[:,1])

    for f,c in zip(phi_array, colors):
        phi = data[np.isclose(data[:,3],f) & np.isclose(data[:,2],10e5),:]
        axes.semilogy(1000/T_array, phi[:,7],ls=ls, color=c, lw=lw, label=f"{pl} phi={f}")



axes.legend()

axes.set_xlabel("1000/T")
axes.set_ylabel("IDT (s)")

plt.show()

