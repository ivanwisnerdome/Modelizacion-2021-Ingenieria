import numpy as np
import matplotlib.pyplot as plt

nx = 50
ny = 50
dX = 1/nx
dY = 1/ny
beta = 1

"""
Condiciones de borde para temperatura en los diferentes lados
"""

TArriba = 100
Tabajo = 20
TIzq = 75
TDer = 50

"""
Condiciones de borde para flujo en los diferentes lados
"""

QArriba = 0
Qabajo = 0
QIzq = 0
QDer = 0

"""
La lista CdC te dice (en orden izquierda, derecha, abajo, arriba) si la
condición que le voy a poner a ese lado es de temperatura (en tal caso
ese indice vale 1), o de flujo (vale 0)
"""

CdC = [1,0,0,1]


def kvalue(i,j,nx):
	K = (i + j*nx)
	return K

def ijvalue(k,nx):
	j = np.floor(k/nx)
	i = k - j*nx
	return i, j

maximo = nx*ny
Malla = np.eye(maximo)
MallaT = np.zeros((ny,nx))
MallaDerx = np.zeros((ny,nx))
MallaDery = np.zeros((ny,nx))
DerX = np.zeros((ny,nx))
DerY = np.zeros((ny,nx))
Temp = np.zeros(maximo)
b = np.zeros(maximo)

for i in range(CdC[0],nx-CdC[1]):
	for j in range(CdC[2], ny-CdC[3]):

		K = kvalue(i,j,nx)
		Malla[K,K] = -2*(1+beta**2)

		if j == 0:
			Malla[K,(K + nx)] = 2*beta**2
			if i == 0:
				Malla[K,(K + 1)]  = 2
			elif i == nx - 1:
				Malla[K,(K - 1)]  = 2
			else:
				Malla[K,(K - 1)]  = 1
				Malla[K,(K + 1)]  = 1
		elif j == ny - 1:
			Malla[K,(K - nx)] = 2*beta**2
			if i == 0:
				Malla[K,(K + 1)]  = 2
			elif i == nx - 1:
				Malla[K,(K - 1)]  = 2
			else:
				Malla[K,(K - 1)]  = 1
				Malla[K,(K + 1)]  = 1
		elif i == 0:
			Malla[K,(K - nx)] = beta**2
			Malla[K,(K + nx)] = beta**2
			Malla[K,(K + 1)]  = 2
		elif i == nx - 1:
			Malla[K,(K - nx)] = beta**2
			Malla[K,(K + nx)] = beta**2
			Malla[K,(K - 1)]  = 2
		else:
			Malla[K,(K + nx)] = beta**2
			Malla[K,(K - nx)] = beta**2
			Malla[K,(K - 1)]  = 1
			Malla[K,(K + 1)]  = 1

for i in range(CdC[2],ny-CdC[3]):
	if CdC[0] == 1:
		b[kvalue(0,i,nx)]= TIzq
	elif CdC[0] == 0:
		b[kvalue(0,i,nx)]= QIzq
	
for i in range(CdC[2],ny-CdC[3]):
	if CdC[1] == 1:
		b[kvalue(nx-1,i,nx)] = TDer
	elif CdC[1] == 0:
		b[kvalue(nx-1,i,nx)] = QDer

for i in range(CdC[0],nx-CdC[1]):
	if CdC[2] == 1:
		b[i] = Tabajo
	elif CdC[2] == 0:
		b[i] = Qabajo

for i in range(maximo - nx + CdC[0], maximo - CdC[1]):
	if CdC[3] == 1:
		b[i] = TArriba
	elif CdC[3] == 0:
		b[i] = QArriba

T = np.linalg.solve(Malla, b)

for i in range(0,ny):
	for j in range(0,nx):
		MallaT[ny-1-i,j] = T[kvalue(j,i,nx)]

MallaT[0,0] = (MallaT[0,1] + MallaT[1,0])/2
MallaT[0,nx-1] = (MallaT[0,nx-2] + MallaT[1,nx-1])/2
MallaT[ny-1,0] = (MallaT[ny-1,1] + MallaT[ny-2,0])/2
MallaT[ny-1,nx-1] = (MallaT[ny-1,nx-2] + MallaT[ny-2,nx-1])/2

#plt.imshow(MallaT, cmap="plasma")
plt.contourf(MallaT, origin = "upper", cmap="plasma")
plt.colorbar()

for j in range(0,nx):
	for i in range(0,ny):
		if i == 0:
			MallaDery[ny-i-1,j] = (MallaT[i,j]-MallaT[i+1,j])/dY
			if j == 0:
				MallaDerx[ny-i-1,j] = (MallaT[i,j+1]-MallaT[i,j])/dX
			elif j == nx - 1:
				MallaDerx[ny-i-1,j] = (MallaT[i,j]-MallaT[i,j-1])/dX
			else:
				MallaDerx[ny-i-1,j] = (MallaT[i,j+1]-MallaT[i,j-1])/(2*dX)
		elif i == ny-1:
			MallaDery[ny-i-1,j] = (MallaT[i-1,j]-MallaT[i,j])/dY
			if j == 0:
				MallaDerx[ny-i-1,j] = (MallaT[i,j+1]-MallaT[i,j])/dX
			elif j == nx - 1:
				MallaDerx[ny-i-1,j] = (MallaT[i,j]-MallaT[i,j-1])/dX
			else:
				MallaDerx[ny-i-1,j] = (MallaT[i,j+1]-MallaT[i,j-1])/(2*dX)
		elif j == 0:
			MallaDerx[ny-i-1,j] = (MallaT[i,j+1]-MallaT[i,j])/dX
			MallaDery[ny-i-1,j] = (MallaT[i-1,j]-MallaT[i+1,j])/(2*dY)
		elif j == nx - 1:
			MallaDerx[ny-i-1,j] = (MallaT[i,j]-MallaT[i,j-1])/dX
			MallaDery[ny-i-1,j] = (MallaT[i-1,j]-MallaT[i+1,j])/(2*dY)
		else:
			MallaDerx[ny-i-1,j] = (MallaT[i,j+1]-MallaT[i,j-1])/(2*dX)
			MallaDery[ny-i-1,j] = (MallaT[i-1,j]-MallaT[i+1,j])/(2*dY)

		DerY[i,j] = i
		DerX[i,j] = j

plt.streamplot(DerX + 0.5, DerY + 0.5, -MallaDerx, -MallaDery, density = 2)
axes = plt.gca()
axes.set_xlim([0.5,nx-0.5])
axes.set_ylim([0.5,ny-0.5])
plt.show()


"""
#plt.contourf(MallaT, origin = "upper",cmap="plasma")

izquierda = []
derecha = []

for i in range(CdC[2],ny-CdC[3]):
	izquierda.append(kvalue(0,i,nx))
	derecha.append(kvalue(nx-1,i,nx))

abajo = list(range(CdC[0],nx-CdC[1]))
arriba = list(range(maximo - nx + CdC[0], maximo - CdC[1]))

if CdC[0] == 1:
	for i in izquierda:
		b[i] = TIzq

if CdC[1] == 1:
	for i in derecha:
		b[i] = TDer

if CdC[2] == 1:
	for i in abajo:
		b[i] = Tabajo

if CdC[3] == 1:
	for i in arriba:
		b[i] = TArriba

"""