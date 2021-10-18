
import sys
import numpy as np
import matplotlib.pyplot as plt

assert(len(sys.argv) > 1), "Must supply a log file."

log_file = sys.argv[1]

data=np.genfromtxt(log_file, delimiter=",")
opt_vals = data[:,-1]
x = np.array(range(len(opt_vals)))+1
ten_best_idxs = opt_vals.argsort()[:10]
x_best = x[ten_best_idxs]
best_vals = opt_vals[ten_best_idxs]

plt.plot(x, opt_vals, 'b-', label="All")
plt.plot(x_best, best_vals, 'ro', label="Ten best")
plt.plot(x_best[0], best_vals[0], 'gs', label="Best")
plt.xlabel("Iteration")
plt.ylabel("Objective Function")
plt.ylim([0, 5*best_vals[0]])
plt.title(f"Minimum score: {best_vals[0]} (iteration {ten_best_idxs[0]})")
plt.legend()
plt.show()

