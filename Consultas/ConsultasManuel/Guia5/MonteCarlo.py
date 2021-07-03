import numpy as np
import matplotlib.pyplot as plt
import matplotlib
np.set_printoptions(precision = 4, linewidth = 150)

nx = 50
N = 5000
Grilla = np.round(np.random.rand(nx,nx))*2-1
#print(Grilla)
Jota = 1
Kboltzmann = 1
Temperatura = 100
AceptadosM = 0
Energia = 0
Magnetizacion = 0
MagnetizacionContador = 0
M = []

def kvalue(i,j,nx):
	K = (i + j*nx)
	return K

def ijvalue(k,nx):
	j = (np.floor(k/nx)).astype(int)
	i = (k - j*nx).astype(int)
	return i, j

for i in range(0,nx):
	for j in range(0,nx):
		Energia += -Jota*(Grilla[i,j]*(Grilla[i-1,j] + Grilla[i,j-1]))
		Magnetizacion += Grilla[i,j]

MagnetizacionContador += Magnetizacion
#print(Magnetizacion)
#print(Energia)

for i in range(N):
	AzarK = (np.round(np.random.rand(1)*(nx*nx-1))).astype(int)
	#print(AzarK)
	Indices = ijvalue(AzarK,nx)
	#print(Indices[0],Indices[1])
	DifEnergia = 2 * Jota * Grilla[Indices[0],Indices[1]]*(
		Grilla[Indices[0]-1,Indices[1]] + Grilla[Indices[0],Indices[1]-1] + \
		Grilla[Indices[0],Indices[1]+1-nx] + Grilla[Indices[0]+1-nx,Indices[1]]
		)
	if DifEnergia < 0:
		Grilla[Indices[0],Indices[1]] = -Grilla[Indices[0],Indices[1]]
		AceptadosM += 1
		Magnetizacion += 2*Grilla[Indices[0],Indices[1]]
	else:
		Probabilidad = np.exp(-DifEnergia/(Kboltzmann*Temperatura))
		NumeroAleat = np.random.rand(1)

		if Probabilidad > NumeroAleat:
			AceptadosM += 1
			Grilla[Indices[0],Indices[1]] = -Grilla[Indices[0],Indices[1]]
			Magnetizacion += 2*Grilla[Indices[0],Indices[1]]
	
	MagnetizacionContador += Magnetizacion
	M.append(MagnetizacionContador/((i+2)*nx**2))

plt.plot(M)
plt.show()
plt.imshow(Grilla, cmap = 'binary')
plt.show()

# -J * (s1*s0 + s1*s2)
#  J * (s1*s0 + s1*s2)
#2*J * (s1*s0 + s1*s2) = s1*2*J*(s0 + s2)