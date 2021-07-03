import numpy as np
import matplotlib.pyplot as plt

intervalo = [0, 30]
f = lambda x: 200*(x/(5+x))*np.exp((-2*x)/300)
# pasos = int(input())

def trapecios (intervalo,f,pasos):

	Espaciado = np.linspace(intervalo[0], intervalo[-1], pasos)
	alfa = f(intervalo[0])

	for i in range(1,len(Espaciado)-1):

		alfa += 2*f(Espaciado[i])

	I = ((intervalo[1] - intervalo[0])/len(Espaciado))*((f(intervalo[-1])+alfa)/2)

	return I

# J = trapecios(intervalo, f, pasos)

# print(J)

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

# I = Simpsons(intervalo, f, pasos)

# print(I)


evaluacion = np.logspace(2,4,20).astype(int)*2
IntegralT = []
IntegralS = []

print(evaluacion)
for i in range(len(evaluacion)):
	IntegralT.append(trapecios(intervalo,f,evaluacion[i]))
	IntegralS.append(Simpsons(intervalo,f,evaluacion[i]))

plt.semilogx(evaluacion,IntegralS,'ro',label="Simpson")
plt.semilogx(evaluacion,IntegralT,'gs',label="Trapecios")
plt.legend()
plt.ylabel("Valor de la integral")
plt.xlabel("Número de intervalos en integral")
plt.show()


# astype cambia el tipo de variable, porque logspace genera rango en decadas, misma cantidad de 
# números entre 100 y 1000 que 1000 y 10000, esto no cae en numeros enteros, pero con astype(int)
# lo obligo que sean numeros enteros (sería np.logspace(1,3,10).astype(int))