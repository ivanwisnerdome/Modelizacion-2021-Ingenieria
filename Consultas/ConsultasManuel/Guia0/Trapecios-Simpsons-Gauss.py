import numpy as np
import matplotlib.pyplot as plt
import copy

intervalo = [0, 30]
f = lambda x: 200*(x/(5+x))*np.exp((-2*x)/300)

def trapecios (intervalo,f,pasos):

	Espaciado = np.linspace(intervalo[0], intervalo[-1], pasos)
	alfa = f(intervalo[0])

	for i in range(1,len(Espaciado)-1):

		alfa += 2*f(Espaciado[i])

	I = ((intervalo[1] - intervalo[0])/len(Espaciado))*((f(intervalo[-1])+alfa)/2)

	return I

def Simpsons (intervalo, f, pasos):
	Espaciado = np.linspace(intervalo[0], intervalo[-1], pasos)
	alfa = f(intervalo[0])

	for i in range(1,len(Espaciado)-1):
		if i % 2 == 0:
			alfa += 4*f(Espaciado[i])
		else:
			alfa += 2*f(Espaciado[i])

	I = ((intervalo[1] - intervalo[0])/(len(Espaciado)-1))*(alfa + f(intervalo[-1]))/3

	return I

def PuntosyPesos (n):
	MPuntos = np.zeros((7,7))
	MPesos = np.zeros((7,7))
	MPuntos[1,0] = -1/np.sqrt(3)
	MPuntos[2,0] = -np.sqrt(0.6)
	MPuntos[3,0] = -0.861136311594953
	MPuntos[4,0] = -0.906179845938664
	MPuntos[5,0] = -0.932469514203152
	MPuntos[6,0] = -0.949107912342759
	MPuntos[2,1] = -0.339981043584856
	MPuntos[3,1] = -0.538469310105683
	MPuntos[4,1] = -0.661209386466265
	MPuntos[5,1] = -0.741531185599394
	MPuntos[3,2] = -0.238619186083197
	MPuntos[4,2] = -0.405845151377397
	MPuntosb = copy.copy(MPuntos)
	MPuntosb = -np.transpose(MPuntosb)
	MPuntos = MPuntosb + MPuntos
	
	MPesos[1,0] = 1
	MPesos[2,0] = 5/9
	MPesos[3,0] = 0.347854845137454
	MPesos[4,0] = 0.236926885056189
	MPesos[5,0] = 0.171324492379170
	MPesos[6,0] = 0.129484966168870
	MPesos[2,1] = 0.652145154862546
	MPesos[3,1] = 0.478628670499366
	MPesos[4,1] = 0.360761573048139
	MPesos[5,1] = 0.279705391489277
	MPesos[3,2] = 0.467913934572691
	MPesos[4,2] = 0.381830050505119

	MPesosb = copy.copy(MPesos)
	MPesosb = np.transpose(MPesosb)
	MPesos = MPesosb + MPesos

	MPesos[0,0] = 2
	MPesos[1,1] = 8/9
	MPesos[2,2] = 0.568888888888889
	MPesos[3,3] = 0.417959183673469

	VPuntos = []
	VPesos = []
	for i in range(0, n):
		VPuntos.append(MPuntos[n-i-1,i])
		VPesos.append(MPesos[n-i-1,i])

	return VPuntos, VPesos


def Gauss (intervalo, f, puntos):

	Pendiente = (intervalo[1]-intervalo[0])/2
	OdO = (intervalo[1]+intervalo[0])/2
	VPuntos, VPesos = PuntosyPesos(puntos)

	I = 0

	for i in range(0, puntos):
		I += Pendiente*f(Pendiente*VPuntos[i]+OdO)*VPesos[i]

	return I

evaluacion = np.logspace(0,1,20).astype(int)*2
IntegralT = []
IntegralS = []
IntegralG = []


for i in range(1,8):
	IntegralG.append(Gauss(intervalo,f,i))

for i in range(len(evaluacion)):
	IntegralT.append(trapecios(intervalo,f,evaluacion[i]))
	IntegralS.append(Simpsons(intervalo,f,evaluacion[i]))

print(IntegralG[-1])
print(IntegralS[-1])
print(IntegralT[-1])

plt.semilogx(evaluacion,IntegralS,'r',label = "Simpson")
plt.semilogx(evaluacion,IntegralT,'g',label = "Trapecios")
plt.semilogx(range(1,8),IntegralG,'b', label = "Gauss")
plt.legend()
plt.ylabel("Valor de la integral")
plt.xlabel("Número de intervalos en integral")
plt.show()

"""
Por cierto, encontre que la razón por que la otra clase veíamos que la de Simpsons convergía muy rápido
en el gráfico era porque el gráfico comenzaba a 200 pasos, no a 1 (el vector evaluación estaba definido
como np.logspace(2,4,20).astype(int)*2)
"""