import numpy as np
import matplotlib.pyplot as plt
import sys
import pylab as py

data= np.loadtxt(sys.argv[1])

t=data[:,0]
x=data[:,1]
y=data[:,2]

# nombre del archivo

nombre="poblaciones_"+ str(data[0,1]) + "_" + str(data[0,2]) +".pdf"

title="Grafica poblaciones X Vs Y (xo ="+ str(data[0,1]) + "_yo =" + str(data[0,2]) + ")"

fig=plt.figure()
plt.plot(x, y)
plt.title(title)
plt.xlabel("$X$")
plt.ylabel("$Y$")
#py.show()

fig.savefig(nombre)
