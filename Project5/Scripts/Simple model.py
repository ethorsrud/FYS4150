import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.special import gamma



files = ["histogram_lambda0.txt","histogram_lambda0.25.txt","histogram_lambda0.5.txt","histogram_lambda0.9.txt"]
legend = [r"$\lambda$ = 0",r"$\lambda$ = 0 Gibbs",r"$\lambda$ = 0.25",r"$\lambda$ = 0.25 P(x)",r"$\lambda$ = 0.5",r"$\lambda$ = 0.5 P(x)",r"$\lambda$ = 0.9",r"$\lambda$ = 0.9 P(x)"]
for file in files:
    text = open(file,"r")

    prob = []

    for line in text:
        prob.append(float(line.split()[0]))

    text.close()

    mean = np.mean(prob)
    N = len(prob)
    beta = 1/mean
    bins = np.linspace(0,500,51)
    hist,be = np.histogram(np.array(prob), bins)
    x = 0.5*(be[1:]+be[:-1])

    omega = beta*np.exp(-beta*x)

    lamb = float(file[16:-4])
    print("Lambda",lamb)
    n = 1+(3*lamb)/(1-lamb)
    an = ((n**n)/gamma(n))
    xarr = np.linspace(np.min(prob),np.max(prob),1000)
    xarr2 = xarr/mean#x/mean
    P = an*xarr2**(n-1)*np.exp(-n*xarr2)

    plt.plot(x,hist/N)

    if file == "histogram_lambda0.txt":
        plt.plot(x,omega*10,'--')
    else:
        plt.plot(xarr,P/10,'--')

    plt.xlabel("m")
    plt.ylabel("Probability")
    plt.title("m0 = 100")
    #plt.semilogy(x,y)
    #plt.xlabel("m")
    #plt.ylabel("Probability")

    #plt.semilogy(x,y)
    #plt.show()
#plt.title("logplot - m0 = 100")
plt.legend(legend)
plt.show()
