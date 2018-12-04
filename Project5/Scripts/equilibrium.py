import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.special import gamma



files = ["mc10.txt","mc50.txt","mc100.txt","mc200.txt"]
legend = ["10 simulations","50 simulations","100 simulations","200 simulations"]
for file in files:
    text = open(file,"r")

    prob = []

    for line in text:
        prob.append(float(line.split()[0]))

    text.close()

    #n, bins, patches = plt.hist(prob,bins = 100,density=True)
    mean = np.mean(prob)
    N = len(prob)
    beta = 1/mean
    bins = np.linspace(0,500,501)
    hist,be = np.histogram(np.array(prob), bins)
    x = 0.5*(be[1:]+be[:-1])

    plt.plot(x,hist/N)


    plt.xlabel("m")
    plt.ylabel("Probability")
    plt.title("m0 = 100")

plt.legend(legend)
plt.show()
