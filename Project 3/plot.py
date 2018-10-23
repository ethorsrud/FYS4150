import numpy as np
import matplotlib.pyplot as plt
import sys

planets = int(sys.argv[1])

files = [str(i)+".txt" for i in range(planets)]

legend = ["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune","Pluto"]

for file in files:
    text = open(file,"r")

    t = []
    x = []
    y = []
    v_x = []
    v_y = []

    for line in text:
        t.append(float(line.split()[0]))
        x.append(float(line.split()[1]))
        y.append(float(line.split()[2]))
        v_x.append(float(line.split()[3]))
        v_y.append(float(line.split()[4]))

    text.close()
    plt.plot(x,y)

plt.axis("equal")

#plt.plot(0,0,'*',color="y",markersize = 15)
plt.xlabel("X")
plt.ylabel("Y")
plt.legend(legend,loc='upper right')
plt.show()
