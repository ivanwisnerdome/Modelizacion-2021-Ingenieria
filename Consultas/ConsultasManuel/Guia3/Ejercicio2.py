import numpy as np
import matplotlib.pyplot as plt
import gmsh
np.set_printoptions(precision = 4, linewidth = 150)

gmsh.initialize()
gmsh.model.add('test')

lc = 0.4
L = 5

p1 = gmsh.model.geo.addPoint(0  ,  L, 0, lc)
p2 = gmsh.model.geo.addPoint(2*L,  0, 0, lc)
p3 = gmsh.model.geo.addPoint(2*L,  L, 0, lc/5)
p4 = gmsh.model.geo.addPoint(0  ,  0, 0, lc)
p5 = gmsh.model.geo.addPoint(0.5,  0, 0, lc/10)
p6 = gmsh.model.geo.addPoint(0  ,0.5, 0, lc/10)
l1 = gmsh.model.geo.addLine(p6,p1)
l2 = gmsh.model.geo.addLine(p2,p5)
l3 = gmsh.model.geo.addLine(p1,p3)
l4 = gmsh.model.geo.addLine(p3,p2)
l5 = gmsh.model.geo.addCircleArc(p5,p4,p6)

C1 = gmsh.model.geo.addCurveLoop([l1, l3, l4, l2, l5])
S1 = gmsh.model.geo.addPlaneSurface([C1])

gmsh.model.geo.synchronize()
EmpotradoX = gmsh.model.addPhysicalGroup(1, [l1])
gmsh.model.setPhysicalName(1, EmpotradoX, 'EmpotradoX')

EmpotradoY = gmsh.model.addPhysicalGroup(1, [l2])
gmsh.model.setPhysicalName(1, EmpotradoY, 'EmpotradoY')

Traccionado = gmsh.model.addPhysicalGroup(1, [l4])
gmsh.model.setPhysicalName(1, Traccionado, 'Traccionado')

Superficie = gmsh.model.addPhysicalGroup(2, [S1])
gmsh.model.setPhysicalName(2, Superficie, 'Superficie')

TrashNodo = gmsh.model.addPhysicalGroup(0, [p4])
gmsh.model.setPhysicalName(0, TrashNodo, 'Trash')

EsquinasTracc = gmsh.model.addPhysicalGroup(0, [p3, p2])
gmsh.model.setPhysicalName(0, EsquinasTracc, 'Esq')

gmsh.model.geo.synchronize()

gmsh.model.mesh.generate(2)
gmsh.model.geo.synchronize()

NodeInfo = gmsh.model.mesh.get_nodes()
NumeroNodos = NodeInfo[0].shape[0]
MatrizNodos = NodeInfo[1].reshape(NumeroNodos, 3)

ElementInfo = gmsh.model.mesh.get_elements()
ETYPES = ElementInfo[0]

ETAGS, ELEMENTS = gmsh.model.mesh.get_elements_by_type(2)
MatrizConectividad = ELEMENTS.reshape([ETAGS.shape[0],3]).astype(int) - 1
NodosEmpotradosX = gmsh.model.mesh.get_nodes_for_physical_group(1,EmpotradoX)[0].astype(int) - 1
NodosEmpotradosY = gmsh.model.mesh.get_nodes_for_physical_group(1,EmpotradoY)[0].astype(int) - 1
NodosTrash = gmsh.model.mesh.get_nodes_for_physical_group(0, TrashNodo)[0].astype(int) - 1

# Las entidades son las cosas asociadas a physical groups, estos solo son geometría

entityTraccionada = gmsh.model.getEntitiesForPhysicalGroup(1, Traccionado)

# Con esto puedo tener los elementos para ese número de entidades

Tgroup, Ttraccionada, Ltraccionada = gmsh.model.mesh.getElements(1, entityTraccionada[0])

# Tengo una matriz de líneas y las longitudes las saco de ahí con la diferencia
# Luego hago un loop para ver cuanta fuerza le asigno a cada nodo
# La fuerza locales tensión espesor longitudes y la fuerza del elemento (creo)
# sobre dos porque reparto entre los dos nodos del elemento

Ltraccionada = np.ravel(Ltraccionada).astype(int) - 1
Length = int(len(Ltraccionada)/2)
elementos = MatrizConectividad.shape[0]
nodos = MatrizNodos.shape[0]
Espesor = 1
u = np.zeros(2*nodos)
F = np.zeros(2*nodos)
FuerzaDistr = 1000
Fuerza = 10000
ParedL = 10
poisson = 0.3
E = 30e6

for i in range(Length):
	Var = [2*i,2*i+1]
	X = MatrizNodos[Ltraccionada[Var],0]
	Y = MatrizNodos[Ltraccionada[Var],1]
	PL = np.sqrt((X[1]-X[0])**2 + (Y[1]-Y[0])**2)
	F[2*Ltraccionada[2*i]] += Fuerza*(PL)/(ParedL*2)
	F[2*Ltraccionada[2*i+1]] += Fuerza*(PL)/(ParedL*2)

S = []

for i in NodosEmpotradosX:
	S.append(2*i)
for i in NodosEmpotradosY:
	S.append(2*i + 1)
for i in NodosTrash:
	S.append(2*i)
	S.append(2*i+1)

S = np.ravel(np.sort(S))

R = [i for i in range(2*nodos) if i not in S]
R = np.ravel(R)

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
D = D*E/(1 - poisson**2)

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

	KLocal = Espesor*np.absolute(A[i])*np.dot(np.transpose(B), np.dot(D,B))

	IndicA = []
	Indic = []

	for m in range(3):
		IndicA.append(np.linspace(2*N[m], (2*N[m] + 1), 2))

	Indic = np.ravel(IndicA).astype(int)

	KGlobal[np.ix_(Indic,Indic)] += KLocal

u[R] = np.linalg.solve(KGlobal[np.ix_(R,R)], F[R] - KGlobal[np.ix_(R, S)].dot(u[S]))
F[S] = KGlobal[S,:].dot(u)

for i in range(elementos):
	N = MatrizConectividad[i,:]

	IndicA = []
	Indic = []
	for j in range(3):
		IndicA.append(np.linspace(2*N[j], (2*N[j] + 1), 2))
	Indic = np.ravel(IndicA).astype(int)

	Tensiones.append(np.dot(D, np.dot(Btot[i], u[Indic])))

TensionesPromM = np.zeros((nodos,2))
Tensiones = np.array(Tensiones)
TensionesProm = np.zeros(nodos)

#MDF-COMMENT no me queda claro esto.
for i in range(elementos):
	N = MatrizConectividad[i,:]
	TensionesPromM[N[0],0] += Tensiones[i,0]
	TensionesPromM[N[1],0] += Tensiones[i,0] #MDF-COMMENT aca no es 
	TensionesPromM[N[2],0] += Tensiones[i,0]
	TensionesPromM[N[0],1] += 1
	TensionesPromM[N[1],1] += 1
	TensionesPromM[N[2],1] += 1

#MDF-COMMENT no me queda claro por  que lo de ariva te funcione.
#MDF-COMMENT for node in range(NumeroNodos):
#MDF-COMMENT     SIGMA_AVE_MAX[node] = SIGMA_MAX[ ETAGS[ (MC == node+1).any(axis=1) ] - ETAGS.min() ].mean()



No2 = [i for i in range(nodos) if i not in NodosTrash]

#MDF-COMMENT no entiendo que querés hacer acá ... 
for i in No2:
	TensionesProm[i] = TensionesPromM[i,0]/TensionesPromM[i,1]

Desplazamientos = np.zeros((nodos,3))
Forces = np.zeros((nodos,3))

for i in No2:
	Desplazamientos[i,0] = u[2*i]
	Desplazamientos[i,1] = u[2*i + 1]

for i in No2:
	Forces[i,0] = F[2*i]
	Forces[i,1] = F[2*i + 1]

desps = gmsh.view.add("Desplazamientos")
Fza = gmsh.view.add("Fuerzas")
ViewTension = gmsh.view.add("Tensión")
TensPromV = gmsh.view.add("Tensión Promedio")

Desps = gmsh.view.addModelData(desps, 0, 'test', 'NodeData', NodeInfo[0], Desplazamientos, numComponents = 3)
Forz = gmsh.view.addModelData(Fza, 0, 'test', 'NodeData', NodeInfo[0], Forces, numComponents = 3)
Tens = gmsh.view.addModelData(ViewTension, 0,  'test', 'ElementData', ETAGS, Tensiones[:,0].reshape(-1,1), numComponents = 1)
TensionesP = gmsh.view.addModelData(TensPromV, 0, 'test', 'NodeData', NodeInfo[0], TensionesProm.reshape(-1,1), numComponents = 1)

gmsh.fltk.run()
gmsh.finalize()
