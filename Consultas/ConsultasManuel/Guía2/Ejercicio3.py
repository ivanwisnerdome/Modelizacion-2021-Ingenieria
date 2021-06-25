import numpy as np
np.set_printoptions(precision = 4, linewidth = 150)
import matplotlib.pyplot as plt

d = lambda x,T,L,A,E: - T*(x**3 - L**3)/(6*A*E)
f = lambda x,T,A: - T*(x**2)/(2*A)

Largo = 1.50
T = -200000
elementos = 10
A = 0.001
E = 210e9

Desplazamientos = np.zeros(elementos + 1)
Fuerzas = np.zeros(elementos + 1)

VectLong = np.linspace(0, Largo, elementos + 1)

KGlob = np.zeros((elementos + 1, elementos + 1))
Kloc = np.eye(2)
Kloc[0,1] = -1
Kloc[1,0] = -1
Tensión = []
Ftot = T*(VectLong[1]**2 - VectLong[0]**2)/2
Fvec = np.zeros(2)
CentroElem = []

for i in range(elementos):

	Fuerzas[i] += Ftot/3 + (VectLong[i + 1] - VectLong[i])*T*VectLong[i]/2
	Fuerzas[i+1] += Ftot*2/3 + (VectLong[i + 1] - VectLong[i])*T*VectLong[i]/2

	K = E*A/(VectLong[i + 1] - VectLong[i])
	KGlob[i : i + 2, i : i + 2] += Kloc*K

Desplazamientos[0:elementos] = np.linalg.solve(KGlob[np.ix_(np.arange(elementos),\
	np.arange(elementos))], Fuerzas[np.arange(elementos)])
Respuesta = np.dot(KGlob[elementos,np.arange(elementos + 1)], Desplazamientos) - Fuerzas[elementos]

for i in range(elementos):
	Tensión.append(E*np.array(
		[-1/(VectLong[i+1] - VectLong[i]), 1/(VectLong[i+1] - VectLong[i])]).dot(Desplazamientos[i:i+2]
								))
	CentroElem.append((VectLong[i + 1] + VectLong[i])/2)

plt.figure(figsize = (12,8))
plt.plot(np.linspace(0, Largo, 100), d(np.linspace(0, Largo, 100),T,Largo,A,E),'r', label = 'Ecuación')
plt.plot(np.linspace(0, Largo, elementos + 1), Desplazamientos, label = 'Aproximación')
plt.legend()
plt.show()

plt.figure(figsize = (12,8))
plt.plot(np.linspace(0, Largo, 100), f(np.linspace(0, Largo, 100),T,A),'r', label = 'Ecuación')
plt.plot(CentroElem, Tensión, 's', label = 'Aproximación')
plt.step(CentroElem, Tensión, 'k', where = 'mid', label = 'Step')
plt.legend()
plt.show()

#Cada elemento tiene tensiones constantes Sigma = E * [-1/L 1/L]*[d1, d2],son constantes porque 
# uso una aproximación tal que la tensión dentro de cada elemento es constante
#
#TENY = np.vstack( ([0.0], self.Sig) )
#        plt.step(self.MN[:, 0], TENY, 'r',
#                where='pre', label='Solución Numérica')