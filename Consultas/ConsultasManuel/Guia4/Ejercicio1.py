import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import gmsh
from scipy import linalg
from scipy.linalg import eigh
np.set_printoptions(precision = 4, linewidth = 150)

Elemm = 2
ElemM = 3

for m in range(Elemm,ElemM):

	elem = m
	nodos = elem + 1
	Densidad = 7850
	E = 210e9
	Long = 1
	L = Long/elem
	Area = 1e-3
	Inercia = 1e-7
	MatMasa = np.zeros((4,4))
	MasaTot = np.zeros((2*nodos,2*nodos))
	KSub = np.zeros((4,4))
	KGlobal = np.zeros((2*nodos,2*nodos))
	MatrizConectividad = np.zeros((elem,2))
	MatrizConectividad[:,0] = np.linspace(0, nodos - 2, elem).astype(int)
	MatrizConectividad[:,1] = np.linspace(1, nodos - 1, elem).astype(int)
	S = [0, 1]
	R = [i for i in range(2*nodos) if i not in S]
	
	for i in range(elem):
		N = MatrizConectividad[i,:]	
	
		#MatMasa[0,0] = 12*Densidad*Area*L/24
		#MatMasa[1,1] = L**2*Densidad*Area*L/24
		#MatMasa[2,2] = 12*Densidad*Area*L/24
		#MatMasa[3,3] = L**2*Densidad*Area*L/24
		
		MatMasa = (Densidad*Area*L/420)*np.array([ 
														[  156,    22*L,    54,    -13*L],
														[ 22*L,  4*L**2,  13*L,  -3*L**2],
														[   54,    13*L,   156,    -22*L],
														[-13*L, -3*L**2, -22*L,   4*L**2],
					  							])
	
		KSub = (E*Inercia/L**3)*np.array([ 
										[ 12,    6*L,  -12,    6*L],
										[6*L, 4*L**2, -6*L, 2*L**2],
										[-12,   -6*L,   12,   -6*L],
										[6*L, 2*L**2, -6*L, 4*L**2],
					  					])
		I = np.linspace(2*N[0], 2*N[1] + 1, 4).astype(int)
		KGlobal[np.ix_(I, I)] += KSub
		MasaTot[np.ix_(I, I)] += MatMasa
	
	Val, Vec = eigh(KGlobal[np.ix_(R,R)], MasaTot[np.ix_(R,R)])
	#print(MasaTot)
	#print(KGlobal)
	Freq = (np.sqrt(Val)/(2*np.pi)).reshape(-1,1)
	#for i in range(0,elem):
	#	plt.plot(nodos,Freq[2*i],'ko')
	#	plt.plot(nodos,Freq[2*i+1],'ro')
	plt.plot(nodos,Freq[0],'ko')
	plt.plot(nodos,Freq[1],'ro')
	plt.plot(nodos,Freq[2],'bo')
	plt.plot(nodos,Freq[3],'go')
	plt.xlabel("# de Nodos", fontsize = 10)
	plt.ylabel("Frecuencia", fontsize = 10)
#print(Vec.reshape(-1,2*(nodos-1)))
#print(Freq)
plt.show()
Utot = np.zeros((Vec.shape[1]+2,Vec.shape[0]))
V = np.arange(2,2*nodos)
Utot[V,:] = Vec

V2 = np.arange(0,nodos)*2
Normaliz = Utot[-2,:]

v = lambda a1,a2,a3,a4,x: a1*x**3 + a2*x**2 + a3*x + a4

for j in range(3):
	Utot[:,j] = Utot[:,j]/Normaliz[j]
	plt.plot(np.linspace(0,1,nodos),Utot[V2,j],'k')
	#for i in range(elem):
	#	A1 = 2/L**3*(Utot[2*i,j] - Utot[2*(i + 1),j]) + 1/L**2*(Utot[2*i + 1,j] + Utot[2*(i + 1) + 1,j])
	#	A2 = - 3/L**2*(Utot[2*i,j] - Utot[2*(i + 1),j]) - 1/L*(2*Utot[2*i + 1,j] + Utot[2*(i + 1) + 1,j])
	#	A3 = Utot[2*i + 1,j]
	#	A4 = Utot[2*i,j]
	#	plt.plot(np.linspace(i*L, (i+1)*L, 25), v(A1, A2, A3, A4, np.linspace(0, L, 25)), 'k')
plt.xlabel("Distancia", fontsize = 10)
plt.ylabel("Amplitud", fontsize = 10)
plt.show()