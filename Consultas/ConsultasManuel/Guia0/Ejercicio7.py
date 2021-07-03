import numpy as np
import matplotlib.pyplot as plt
import copy

f = lambda x, y: 0.5*y
pasos = 10
intervalo = [0, 4]
y0 = 2


def RK(f, pasos, y0, intervalo):
	y = [y0]
	x = np.linspace(intervalo[0], intervalo[1], pasos +1)
	k = [0, 0, 0, 0]
	Espaciado = (intervalo[1] - intervalo[0])/pasos

	for i in range(0,pasos):
		k[0] = f(x[i], y[i])
		k[1] = f(x[i] + Espaciado/2 , y[i] + k[0]*Espaciado/2)
		k[2] = f(x[i] + Espaciado/2 , y[i] + k[1]*Espaciado/2)
		k[3] = f(x[i] + Espaciado, y[i] +k[2]*Espaciado)
		y.append(y[i] + 1/6*(k[0] + 2*k[1] + 2*k[2] + k[3])*Espaciado)

	return y

def Euler(f, pasos, y0, intervalo):
	y = [y0]
	x = np.linspace(intervalo[0], intervalo[1], pasos +1)
	Espaciado = (intervalo[1] - intervalo[0])/pasos

	for i in range(0,pasos):
		y.append(y[i] + f(x[i], y[i])*Espaciado)

	return y

def Heun(f, pasos, y0, intervalo):
	y = [y0]
	x = np.linspace(intervalo[0], intervalo[1], pasos +1)
	Espaciado = (intervalo[1] - intervalo[0])/pasos

	for i in range(0,pasos):
		y.append(y[i] + (f(x[i], y[i]) + f(x[i] + Espaciado, y[i] + f(x[i], y[i])*Espaciado))*Espaciado/2)

	return y

"""

g = lambda x: 4/1.3*(np.exp(0.8*x)-np.exp(-0.5*x))+2*np.exp(-0.5*x)

plt.plot(x,y,'r')

X = np.linspace(intervalo[0], intervalo[1], 200)

plt.plot(X,g(X),'b')

plt.show()

"""

ErrR = []
ErrE = []
ErrH = []
maximo = 500 #MDF-COMMENT 50 son pocos puntos

#MDF-COMMENT el tema del rango lineal es que tenes mucha información 'repetida', es mejor hacerlo en escala logarítmica:
nintervalos = np.logspace(1,np.log10(maximo), 10).astype(int)
for i in nintervalos: #range(1,maximo):
	yR = RK(f, i, y0, intervalo)
	yE = Euler(f, i, y0, intervalo)
	yH = Heun(f, i, y0, intervalo)
	ErrR.append(abs(75.3389626-yR[-1])/75.3389626)
	ErrE.append(abs(75.3389626-yE[-1])/75.3389626)
	ErrH.append(abs(75.3389626-yH[-1])/75.3389626)
#MDF-COMMENT yo lo pondría con escala logarítmica
#plt.plot(range(1,maximo),ErrR,'r', label = "Runge-Kutta")
#plt.plot(range(1,maximo),ErrE,'g', label = "Euler")
#plt.plot(range(1,maximo),ErrH,'b', label = "Heun")
plt.loglog(nintervalos,ErrR,'r', label = "Runge-Kutta")
plt.loglog(nintervalos,ErrE,'g', label = "Euler")
plt.loglog(nintervalos,ErrH,'b', label = "Heun")
#MDF-COMMENT para entender lo que graficás, le tenes que poner nombre a los ejes !
plt.ylabel('Error')
plt.xlabel('numero de intervalos')

plt.legend()
#MDF-COMMENTplt.show()
plt.savefig('ResultadoRaro.pdf')
