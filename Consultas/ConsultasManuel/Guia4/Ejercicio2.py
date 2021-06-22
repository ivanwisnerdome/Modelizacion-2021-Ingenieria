import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import gmsh
from scipy import linalg
from scipy.linalg import eigh
np.set_printoptions(precision = 4, linewidth = 150)

elem = 4
nodos = elem + 1
Longitud = 0.1
L = Longitud/elem
Diametro = 0.012
Temperatura0 = 30
TemperaturaIzq = 80
Conductividad = 200
CalorEsp = 900
Densidad = 2700
Intervalot = 0.1

MatCap = np.zeros((2,2))
CapTot = np.zeros((nodos,nodos))
KSub = np.zeros((2,2))
KGlobal = np.zeros((nodos,nodos))
MatrizConectividad = np.zeros((elem,2))
MatrizConectividad[:,0] = np.linspace(0, nodos - 2, elem).astype(int)
MatrizConectividad[:,1] = np.linspace(1, nodos - 1, elem).astype(int)
S = [0,nodos-1]
R = [i for i in range(nodos) if i not in S]

for i in range(elem):
	N = MatrizConectividad[i,:]	

	MatCap[0,0] = 2
	MatCap[1,1] = 2
	MatCap[0,1] = 1
	MatCap[1,0] = 1
	MatCap = CalorEsp*Densidad*L/6*MatCap
	
	KSub = (Conductividad/L)*np.array([
										[1,-1],
										[-1,1]
									])

	I = np.linspace(N[0], N[1], 2).astype(int)
	KGlobal[np.ix_(I, I)] += KSub
	CapTot[np.ix_(I, I)] += MatCap

T = np.zeros(nodos)
T[0] = TemperaturaIzq
T[R] = Temperatura0
T[-1] = Temperatura0
#plt.plot(np.linspace(0, 1, nodos),T)
Contador = 0
Flujo = np.zeros(nodos)

while Contador < 10000:

	Contador += 1
	Tb = T

	T[R] = np.linalg.solve(CapTot[np.ix_(R, R)], np.dot(CapTot[np.ix_(R, R)], T[R]) -\
	 (np.dot(KGlobal[np.ix_(R, R)], T[R]) + np.dot(KGlobal[np.ix_(R, S)], T[S]))*Intervalot)

	Flujo = np.dot(CapTot, (T - Tb)/Intervalot) + np.dot(KGlobal,Tb)

	Error = abs((abs(Flujo[0])-abs(Flujo[-1])))
	#plt.plot(Contador, Error, 'ko')
	plt.plot(Contador, Flujo[0], 'go')
	plt.plot(Contador, -Flujo[-1], 'ro')

	if Error < 1e-4:
		print(Contador)
		break

#plt.plot(np.linspace(0, 1, nodos),T)
plt.show()
