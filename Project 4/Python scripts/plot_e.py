import numpy as np
import matplotlib.pyplot as plt
import sys



files = ["model_1.txt","model_1_up.txt","model_24.txt","model_24_up.txt"]

legend = ["T = 1, random","T = 1, all up","T = 2.4, random","T = 2.4, all up"]

for file in files:
    text = open(file,"r")

    steps = []
    e = []

    for line in text:
        steps.append(float(line.split()[0]))
        e.append(float(line.split()[1]))

    text.close()
    plt.plot(steps,e)

plt.xlabel("MC-steps")
plt.ylabel("Mean energy <E>")
plt.legend(legend,loc='upper right')
plt.show()
