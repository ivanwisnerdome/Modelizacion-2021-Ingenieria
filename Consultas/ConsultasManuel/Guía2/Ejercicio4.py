import numpy as np
np.set_printoptions(precision = 4, linewidth = 150)
import matplotlib.pyplot as plt

u = np.zeros(8)
F = np.zeros(8)
L = 3
elementos = 3
TipoDElem = [1,1,0]

F[4] = -5e4
E = 210e9
I = 2e-4

S = [0, 1, 2, 6, 7]
R = [3, 4, 5]

k = [E*I/(L**3), E*I/(L**3), 200000]

MatrizNodos = np.array([0,3,6,9])
MatrizConectividad = np.array([[0,1],[1,2],[2,3]])
KGlobal = np.zeros((8, 8))

KSub = np.zeros((4,4))

for i in range(elementos):

	N = MatrizConectividad[i,:]

	if TipoDElem[i] == 1:
		KSub = k[i]*np.array([ 
						[ 12,    6*L,  -12,    6*L],
						[6*L, 4*L**2, -6*L, 2*L**2],
						[-12,   -6*L,   12,   -6*L],
						[6*L, 2*L**2, -6*L, 4*L**2],
			  		])

	elif TipoDElem[i] == 0:
		KSub = k[i]*np.array([
						[ 1, 0,-1, 0],
	        			[ 0, 0, 0, 0],
	        			[-1, 0, 1, 0],
	        			[ 0, 0, 0, 0]
	        		])
		
	I = np.linspace(2*N[0], 2*N[1] + 1, 4).astype(int)
	KGlobal[np.ix_(I, I)] += KSub

print(KGlobal/KGlobal.max())

u[R] = np.linalg.solve(KGlobal[np.ix_(R, R)], F[R] - KGlobal[np.ix_(R, S)].dot(u[S]))
F[S] = KGlobal[S,:].dot(u)

print(F)
print(u)

v = lambda a1,a2,a3,a4,x: a1*x**3 + a2*x**2 + a3*x + a4

for i in range(2):
	A1 = 2/L**3*(u[2*i] - u[2*(i + 1)]) + 1/L**2*(u[2*i + 1] + u[2*(i + 1) + 1])
	A2 = - 3/L**2*(u[2*i] - u[2*(i + 1)]) - 1/L*(2*u[2*i + 1] + u[2*(i + 1) + 1])
	A3 = u[2*i + 1]
	A4 = u[2*i]
	plt.plot(np.linspace(0 + 3*i, 3 + 3*i, 30), v(A1, A2, A3, A4, np.linspace(0, L, 30)),'k')

plt.show()