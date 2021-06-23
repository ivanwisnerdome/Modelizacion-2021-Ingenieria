import numpy as np
import matplotlib.pyplot as plt
import gmsh
np.set_printoptions(precision = 4, linewidth = 150)

#MatrizNodosB = np.array([[0, 0], [0, 10], [20, 10], [20, 0], [10, 5]])
#MatrizConectividadB = np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4], [0, 3, 4]])

#MatrizNodos = np.array([[0, 0], [0, 10], [20, 10], [20, 0]])
#MatrizConectividad = np.array([[0, 1, 2], [2, 3, 0]])

MatrizNodos = np.loadtxt("MatrizNodos2D.dat")
MatrizConectividad = np.loadtxt("MatrizConectividad2D.dat", dtype = int) - 1
TraccionadoX = np.loadtxt("TraccionadosX.dat", dtype = int) - 1
Empotrado = np.loadtxt("Empotrados.dat", dtype = int) - 1
EsquinasTracc = np.loadtxt("EsqTracc.dat", dtype = int) - 1

elementos = MatrizConectividad.shape[0]
nodos = MatrizNodos.shape[0]
Espesor = 1
u = np.zeros(2*nodos)
F = np.zeros(2*nodos)
FuerzaDistr = 1000
Fuerza = 10000
poisson = 0.3
E = 30e6

#print(MatrizConectividad)
#print(MatrizNodos)

S = []

for i in Empotrado:
	S.append([2*i, 2*i + 1])

S = np.ravel(S)

R = [i for i in range(2*nodos) if i not in S]
R = np.ravel(R)


#MDF-COMMENT acordate que en realidad hay que medir las longitudes !
for i in TraccionadoX:
	F[2*i] = Fuerza/(len(TraccionadoX)-1)
for i in EsquinasTracc:
	F[2*i] -= Fuerza/(2*(len(TraccionadoX)-1))

print(F)

KLocal = np.zeros((6, 6))
KGlobal = np.zeros((2*nodos, 2*nodos))

A = []
MatAr = np.zeros((3,3))
MatAr[:,0] = 1
Tensiones = []
Btot = []

D = np.zeros((3,3))
D[0,0] = 1
D[1,1] = 1
D[2,2] = 0.5*(1 - poisson)
D[0,1] = poisson
D[1,0] = poisson
D = D*E/(1-poisson**2)

#MDF-COMMENT ya que estas guardando arhivos, podes saltearte todos los cáclulos si ya los guardaste o 
#MDF-COMMENT podes forzarlo a que los haga:

import os

RESULTADO_DESPLAZAMIENTO = 'Desplazamientos.dat'
RESULTADO_TENSIONES = 'Tensiones.dat'
KGLOBAL = 'Kglobal.dat'
Force = False

if os.path.exists(KGLOBAL):
    KGlobal = np.loadtxt(KGLOBAL)
else:
    for i in range(elementos):
            B = np.zeros((3,6))
            beta = np.zeros(3)
            gamma = np.zeros(3)
            
            N = MatrizConectividad[i,:]
            X = MatrizNodos[N,0]
            Y = MatrizNodos[N,1]

            MatAr[:,1] = X
            MatAr[:,2] = Y

            A.append(np.linalg.det(MatAr)/2)

            beta[0] = Y[1] - Y[2]
            beta[1] = Y[2] - Y[0]
            beta[2] = Y[0] - Y[1]
            gamma[0] = X[2] - X[1]
            gamma[1] = X[0] - X[2]
            gamma[2] = X[1] - X[0]
            
            for j in range(3):
                    B[0,2*j] = beta[j]
                    B[1,2*j + 1] = gamma[j]
                    B[2,2*j] = gamma[j]
                    B[2,2*j + 1] = beta[j]

            B = B/(2*A[i])
            Btot.append(B)

            KLocal = Espesor*np.absolute(A[i])*np.dot(np.transpose(B),np.dot(D,B))

            IndicA = []
            Indic = []

            for m in range(3):
                    IndicA.append(np.linspace(2*N[m], (2*N[m] + 1), 2))

##         clap      clap    clap    clap   ## 
            Indic = np.ravel(IndicA).astype(int)

##         clap      clap    clap    clap   ## 
            KGlobal[np.ix_(Indic,Indic)] += KLocal
    np.savetxt(KGLOBAL, KGlobal, fmt='%1.5f')

#print(KGlobal/(KGlobal.max()))
#print(KGlobal*(0.91/375000))

if os.path.exists(RESULTADO_DESPLAZAMIENTO) and not Force:
    Desplazamientos = np.loadtxt(RESULTADO_DESPLAZAMIENTO) #np.zeros((nodos,3))
    u = Desplazamientos[:,:2].reshape(2*nodos, 1, order='C').ravel()
else:
    u[R] = np.linalg.solve(KGlobal[np.ix_(R,R)], F[R] - KGlobal[np.ix_(R, S)].dot(u[S]))
    Desplazamientos = np.zeros((nodos,3))
    for i in range(nodos):
            Desplazamientos[i,0] = u[2*i]
            Desplazamientos[i,1] = u[2*i + 1]
    np.savetxt('Desplazamientos.dat', Desplazamientos, fmt='%1.5f')

F[S] = KGlobal[S,:].dot(u)

#MDF-COMMENT calculamos las fuerzas:
Fx = F[::2]
Fy = F[1::2]


#print(F)
#print(u)

if os.path.exists(RESULTADO_TENSIONES):
    Tensiones = np.loadtxt(RESULTADO_TENSIONES)
else:
    for i in range(elementos):
            N = MatrizConectividad[i,:]

            IndicA = []
            Indic = []
            for j in range(3):
                    IndicA.append(np.linspace(2*N[j], (2*N[j] + 1), 2))
            Indic = np.ravel(IndicA).astype(int)

            Tensiones.append(np.dot(D, np.dot(Btot[i], u[Indic])))
    np.savetxt(RESULTADO_TENSIONES, np.array( Tensiones ), fmt ='%1.5f')

print(Tensiones)
#print(Tensiones)
#MDF-COMMENT  Desplazamientos = np.zeros((nodos,3))
#MDF-COMMENT  for i in range(nodos):
#MDF-COMMENT  	Desplazamientos[i,0] = u[2*i]
#MDF-COMMENT  	Desplazamientos[i,1] = u[2*i + 1]

#print(np.array(Tensiones))

#MDF-COMMENT np.savetxt('Desplazamientos.dat', Desplazamientos, fmt='%1.5f')
#MDF-COMMENT np.savetxt('Tensiones.dat', np.array(Tensiones), fmt='%1.5f')
#MDF-COMMENT ya que estamos:

import pdb
fig, ax = plt.subplots(1,1, figsize=(16, 8))
# ax.triplot(MatrizNodos[:,0], MatrizNodos[:,1], MatrizConectividad)
tcolors = ax.tripcolor(MatrizNodos[:,0], MatrizNodos[:,1], MatrizConectividad, facecolors = Tensiones[:,0], cmap = 'rainbow')
vectors = ax.quiver(
        MatrizNodos[:,0], 
        MatrizNodos[:,1], 
        Desplazamientos[:,0], 
        Desplazamientos[:,1], 
        np.sqrt(Desplazamientos[:,0]**2 + Desplazamientos[:,1]**2),
        scale=2e-4,  
        units='x',
        cmap = 'winter'
        )
fzas = ax.quiver(
        MatrizNodos[:,0], 
        MatrizNodos[:,1], 
        Fx, Fy, 
        np.sqrt(Fx**2 + Fy**2),
        scale=2000, 
        units='x', 
        cmap = 'Wistia')
cbar0 = plt.colorbar(fzas)
cbar0.set_label('Fuerzas (N)', fontsize=24)
cbar1 = plt.colorbar(vectors)
cbar1.set_label('Desplazamientos (m)', fontsize=24)
cbar2 = plt.colorbar(tcolors)
cbar2.set_label('$\sigma _x$ (Pa)', fontsize=24)
ax.set_xlim([MatrizNodos.min()-3, MatrizNodos.max()+5])
fig.tight_layout()
fig.show()
plt.savefig('Resultado.pdf')

