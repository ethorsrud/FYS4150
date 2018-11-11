import numpy as np
import matplotlib.pyplot as plt
import sys



files = ["Task_b.txt"]


for file in files:
    text = open(file,"r")

    e = []
    e_a = []
    m = []
    m_a = []
    cv = []
    cv_a = []
    x = []
    x_a = []

    for line in text:
        e.append(float(line.split()[0]))
        e_a.append(float(line.split()[1]))
        m.append(float(line.split()[2]))
        m_a.append(float(line.split()[3]))
        cv.append(float(line.split()[4]))
        cv_a.append(float(line.split()[5]))
        x.append(float(line.split()[6]))
        x_a.append(float(line.split()[7]))

    text.close()

steps = np.linspace(100,1000000,10000)
legend = ["Calculated","Analytical"]

plt.plot(steps,e)
plt.plot(steps,e_a)
plt.xlabel("MC-steps")
plt.ylabel("Mean energy <E>")
plt.legend(legend,loc='upper right')

plt.figure()
plt.plot(steps,m)
plt.plot(steps,m_a)
plt.xlabel("MC-steps")
plt.ylabel("Mean absolute magnetization <|M|>")
plt.legend(legend,loc='upper right')

plt.figure()
plt.plot(steps,cv)
plt.plot(steps,cv_a)
plt.xlabel("MC-steps")
plt.ylabel("Heat capacity C_v")
plt.legend(legend,loc='upper right')

plt.figure()
plt.plot(steps,x)
plt.plot(steps,x_a)
plt.xlabel("MC-steps")
plt.ylabel("Magnetic Susceptibility X")
plt.legend(legend,loc='upper right')

plt.show()
