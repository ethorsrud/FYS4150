import numpy as np
import matplotlib.pyplot as plt
import sys



files = ["L40.txt","L60.txt","L80.txt","L100.txt"]

legend = ["L = 40","L = 60","L = 80","L = 100"]

for file in files:
    text = open(file,"r")

    T = []
    e = []
    M = []
    cv = []
    x = []

    for line in text:
        T.append(float(line.split()[0]))
        e.append(float(line.split()[1]))
        M.append(float(line.split()[2]))
        cv.append(float(line.split()[3]))
        x.append(float(line.split()[4]))

    text.close()
    plt.plot(T,x)

plt.legend(legend,loc='upper right')
plt.xlabel("Temperature T")
plt.ylabel("Magnetic susceptibility X")
plt.show()
