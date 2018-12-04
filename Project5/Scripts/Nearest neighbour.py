import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.special import gamma



#files = ["histogram_lambda05_alpha05.txt","histogram_lambda05_alpha10.txt","histogram_lambda05_alpha15.txt","histogram_lambda05_alpha20.txt"]
#files = ["histogram_lambda00_alpha05.txt","histogram_lambda00_alpha10.txt","histogram_lambda00_alpha15.txt","histogram_lambda00_alpha20.txt"]
files = ["histogram_lambda00_alpha05_n1000.txt","histogram_lambda00_alpha10_n1000.txt","histogram_lambda00_alpha15_n1000.txt","histogram_lambda00_alpha20_n1000.txt"]

legend = [r"$\alpha$ = 0.5",r"$\alpha$ = 1.0",r"$\alpha$ = 1.5",r"$\alpha$ = 2.0"]
for file in files:
    text = open(file,"r")

    prob = []

    for line in text:
        prob.append(float(line.split()[0]))

    text.close()
    alpha = float(file[24]+"."+file[25])
    print(alpha)
    mean = np.mean(prob)
    N = len(prob)
    beta = 1/mean
    bins = np.linspace(0,5000,5000)
    hist,be = np.histogram(np.array(prob), bins)

    x = 0.5*(be[1:]+be[:-1])
    m = np.linspace(np.min(x),np.max(x),5000)

    #plt.plot(x,hist/N,markersize=0.2,alpha = 1)
    #plt.plot(x,hist/N)
    plt.plot(m,m**(-1-alpha),"--")

    plt.xlabel("m")
    plt.ylabel("Probability")
    plt.title(r"$m_0$ = 100, $\lambda$ = 0.0, N = 1000")

plt.legend(legend)
plt.show()
