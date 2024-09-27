import numpy as np
np.set_printoptions(precision = 4, linewidth = 150)
import matplotlib.pyplot as plt

F, u, vinculos = np.loadtxt("DatosFuV.dat", delimiter=',	', unpack = True)
E, A = np.loadtxt("DatosEA.dat", delimiter=',	', unpack = True)
MatrizNodos = np.loadtxt("MatrizNodos.dat", delimiter=',	')
MatrizConectividad = np.loadtxt("MatrizConectividad.dat", delimiter=',	')
ElemNod = 2

nodos, dimensiones = np.shape(MatrizNodos)
elementos = np.shape(MatrizConectividad)[0]
Tensiones = np.zeros(elementos)

MatrizConectividad = MatrizConectividad.astype(int)

Klocal = np.zeros((dimensiones, dimensiones))
KlocalB = np.zeros((dimensiones*2, dimensiones*2))
KGlobal = np.zeros((dimensiones*nodos, dimensiones*nodos))
k = []
#LargoV = np.zeros((elementos,2))
angulo = []

for i in range(elementos):
	x1, y1 = MatrizNodos[MatrizConectividad[i,0],:]
	x2, y2 = MatrizNodos[MatrizConectividad[i,1],:]
	N = MatrizConectividad[i,:]
	Largo = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
	k.append(E[i]*A[i]/Largo)
	#k = E[i]*A[i]/Largo
	angulo.append(np.arctan2(y2-y1, x2-x1))
	sen = np.sin(angulo[i])
	cos = np.cos(angulo[i])

	#LargoV[i,0] = Largo

	Klocal[0,0] = cos**2
	Klocal[1,1] = sen**2
	Klocal[0,1] = cos*sen
	Klocal[1,0] = cos*sen
	Klocal = Klocal*k[i]

	KlocalB[0:dimensiones,0:dimensiones] = Klocal
	KlocalB[0:dimensiones,dimensiones:] = -Klocal
	KlocalB[dimensiones:,0:dimensiones] = -Klocal
	KlocalB[dimensiones: ,dimensiones: ] = Klocal

	IndicA = []
	Indic = []

	for j in range(ElemNod):
		IndicA.append(np.linspace(dimensiones*N[j], (dimensiones*N[j] + dimensiones - 1), dimensiones))
	
	Indic = np.ravel(IndicA).astype(int)

	KGlobal[np.ix_(Indic,Indic)] += KlocalB

	#KGlobal[2*N1:(2*N1 + 2), 2*N1:(2*N1 + 2)] += Klocal
	#KGlobal[2*N2:(2*N2 + 2), 2*N2:(2*N2 + 2)] += Klocal
	#KGlobal[2*N1:(2*N1 + 2), 2*N2:(2*N2 + 2)] += -Klocal
	#KGlobal[2*N2:(2*N2 + 2), 2*N1:(2*N1 + 2)] += -Klocal

	plt.plot([x1, x2], [y1, y2], 'k', '-')
	plt.text((x2 + x1)/2, (y2 + y1)/2, i, fontsize = 10, color = 'red')


SubF = []
R = [] #Los nodos donde tengo fuerzas dadas
S = [] #Nodos con desplazamientos conocidos

for i in range(dimensiones*nodos):
	if vinculos[i] == 1:
		SubF.append(F[i])
		R.append(i)
	else:
		S.append(i)

#u[R] = np.linalg.solve(KGlobal[np.ix_(R,R)], F[R])
u[R] = np.linalg.solve(KGlobal[np.ix_(R,R)], F[R] - KGlobal[np.ix_(R, S)].dot(u[S]))
#La última parte permite ingresar condicinoes de contorno para las 'u' diferentes a cero
F[S] = KGlobal[S,:].dot(u)

U = []
Floc = np.zeros(4)
Kext = np.eye(4)

for i in range(elementos):
	x1, y1 = MatrizNodos[MatrizConectividad[i,0],:]
	x2, y2 = MatrizNodos[MatrizConectividad[i,1],:]
	N1, N2 = MatrizConectividad[i,:]

	x1 += u[2*N1]
	x2 += u[2*N2]
	y1 += u[2*N1+1]
	y2 += u[2*N2+1]

	U = [u[2*N1], u[2*N1+1], u[2*N2], u[2*N2+1]]

	#Largo = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
	#LargoV[i,1] = Largo

	#Tensiones[i] = (np.dot([F[2*N2] - F[2*N1], F[2*N2 + 1] - F[2*N1 + 1]], [(x2-x1),(y2-y1)]))/A[i]
	x1 += u[2*N1]*149
	x2 += u[2*N2]*149
	y1 += u[2*N1+1]*149
	y2 += u[2*N2+1]*149

	sen = np.sin(angulo[i])
	cos = np.cos(angulo[i])

	L = np.arange(4)

	Kext[L,L] = cos
	Kext[0,1] = sen
	Kext[1,0] = -sen
	Kext[2,3] = sen
	Kext[3,2] = -sen

	Floc = np.dot(k[i]*Kext, U)
	Tensiones[i] = (Floc[2] - Floc[0])/A[i]
	
	plt.plot([x1, x2], [y1, y2], 'b', '-')



for i in range(nodos):
	plt.text(MatrizNodos[i,0],MatrizNodos[i,1], i, fontsize = 14)

print(Tensiones)
print(u)
print(F)
#print(F)
#print(u)
#print(LargoV)
plt.legend(["Antes", "Después"])
plt.show()
#plt.savefig('Puente_Manuel.pdf')
plt.close()