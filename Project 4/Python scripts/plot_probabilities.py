import numpy as np
import matplotlib.pyplot as plt
import sys



files = ["prob_24.txt"]

legend = ["T = 2.4"]
for file in files:
    text = open(file,"r")

    prob = []

    for line in text:
        prob.append(float(line.split()[0]))

    text.close()

    plt.hist(prob,56,weights=(np.ones_like(prob)/float(len(prob))))

plt.xlabel("Energies E")
plt.ylabel("Probability P(E)")
plt.legend(legend,loc='upper right')
plt.show()
