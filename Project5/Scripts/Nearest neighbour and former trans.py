import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.special import gamma



#files = ["N1000_lamb00_alph10_gam00.txt","N1000_lamb00_alph10_gam10.txt","N1000_lamb00_alph10_gam20.txt","N1000_lamb00_alph10_gam30.txt","N1000_lamb00_alph10_gam40.txt"]
#files = ["N1000_lamb00_alph20_gam00.txt","N1000_lamb00_alph20_gam10.txt","N1000_lamb00_alph20_gam20.txt","N1000_lamb00_alph20_gam30.txt","N1000_lamb00_alph20_gam40.txt"]
files = ["N1000_lamb05_alph10_gam00.txt","N1000_lamb05_alph10_gam10.txt","N1000_lamb05_alph10_gam20.txt","N1000_lamb05_alph10_gam30.txt","N1000_lamb05_alph10_gam40.txt"]
#files = ["N1000_lamb05_alph20_gam00.txt","N1000_lamb05_alph20_gam10.txt","N1000_lamb05_alph20_gam20.txt","N1000_lamb05_alph20_gam30.txt","N1000_lamb05_alph20_gam40.txt"]

legend = [r"$\alpha$ = 1.0, $\gamma$ = 0.0",r"$\alpha$ = 1.0, $\gamma$ = 1.0",r"$\alpha$ = 1.0, $\gamma$ = 2.0",r"$\alpha$ = 1.0, $\gamma$ = 3.0",r"$\alpha$ = 1.0, $\gamma$ = 4.0"]
#legend = [r"$\alpha$ = 2.0, $\gamma$ = 0.0",r"$\alpha$ = 2.0, $\gamma$ = 1.0",r"$\alpha$ = 2.0, $\gamma$ = 2.0",r"$\alpha$ = 2.0, $\gamma$ = 3.0",r"$\alpha$ = 2.0, $\gamma$ = 4.0"]

for file in files:
    text = open(file,"r")

    prob = []

    for line in text:
        prob.append(float(line.split()[0]))

    text.close()

    mean = np.mean(prob)
    N = len(prob)
    beta = 1/mean
    bins = np.linspace(0,500,501)
    hist,be = np.histogram(np.array(prob), bins)
    m = np.linspace(0,500,1000)
    x = 0.5*(be[1:]+be[:-1])

    plt.plot(x,hist/N)
    plt.xlabel("m")
    plt.ylabel("Probability")
    plt.title(r"$m_0$ = 100, $\lambda$ = 0.5, N = 1000")

plt.legend(legend)
plt.show()
