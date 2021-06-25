import numpy as np
import matplotlib.pyplot as plt
import gmsh
np.set_printoptions(precision = 4, linewidth = 150)


gmsh.initialize()
gmsh.model.add('test')

lc = 8
L = 10

#Primero puntos
p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
p2 = gmsh.model.geo.addPoint(2*L, 0, 0, lc)
p3 = gmsh.model.geo.addPoint(2*L, L, 0, lc)
p4 = gmsh.model.geo.addPoint(0, L, 0, lc)

#Luego líneas
l1 = gmsh.model.geo.addLine(p1,p2)
l2 = gmsh.model.geo.addLine(p2,p3)
l3 = gmsh.model.geo.addLine(p3,p4)
l4 = gmsh.model.geo.addLine(p4,p1)

#Luego bordes
C1 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
S1 = gmsh.model.geo.addPlaneSurface([C1])

gmsh.model.geo.synchronize()
#Los physical groups son categorías, tengo que darle la
#dimensión a cada grupo y por qué está compuesto
Empotrado = gmsh.model.addPhysicalGroup(1, [l4])
gmsh.model.setPhysicalName(1, Empotrado, 'Empotrado')
Traccionado = gmsh.model.addPhysicalGroup(1, [l2])
gmsh.model.setPhysicalName(1, Traccionado, 'Traccionado')
Superficie = gmsh.model.addPhysicalGroup(2, [S1])
gmsh.model.setPhysicalName(2, Superficie, 'Superficie')
EsquinasTracc = gmsh.model.addPhysicalGroup(0, [p2, p3])
gmsh.model.setPhysicalName(0, EsquinasTracc, 'Esq')

gmsh.model.geo.synchronize()

#Genero el mallado, y quiero que sea 2D
gmsh.model.mesh.generate(2)
#Veo como va quedando la cosa
gmsh.model.geo.synchronize()
#gmsh.fltk.run()

#Obtengo info
NodeInfo = gmsh.model.mesh.get_nodes()
NumeroNodos = NodeInfo[0].shape[0]

MatrizNodos = NodeInfo[1].reshape(NumeroNodos, 3)

#plt.rc('figure', figsize=(10,10))
#plt.plot(MatrizNodos[:,0], MatrizNodos[:,1], 'ok')
#plt.show()

ElementInfo = gmsh.model.mesh.get_elements()
ETYPES = ElementInfo[0]

ETAGS, ELEMENTS = gmsh.model.mesh.get_elements_by_type(2)
MatrizConectividad = ELEMENTS.reshape([ETAGS.shape[0],3])

NodosEmpotrados = gmsh.model.mesh.get_nodes_for_physical_group(1,Empotrado)
NodosTraccionados = gmsh.model.mesh.get_nodes_for_physical_group(1, Traccionado)
NodosEsq = gmsh.model.mesh.get_nodes_for_physical_group(0, EsquinasTracc)

np.savetxt('MatrizNodos2D.dat', MatrizNodos, fmt='%1.3f')
np.savetxt('MatrizConectividad2D.dat', MatrizConectividad, fmt = '%d')
np.savetxt('TraccionadosX.dat', NodosTraccionados[0], fmt = '%d')
np.savetxt('Empotrados.dat', NodosEmpotrados[0], fmt = '%d')
np.savetxt('EsqTracc.dat', NodosEsq[0], fmt = '%d')

DESP = np.loadtxt("Desplazamientos.dat")
Tensiones = np.loadtxt("Tensiones.dat")

#Desplazamiento en función de la coordenada x
#plt.plot(MatrizNodos[:,0], DESP[:,0])
#plt.show()

#Para agregar al modelo tengo que inicial una visualización
desps = gmsh.view.add("desp")
ViewTension = gmsh.view.add("Tension")

#Guardar los datos en la visualización, test es el nombre del modelo, NodeData es para decir que es info de nodos
#porque la info de desplazamiento está asociado a los nodos del mallado, lo de NodeInfo son los números de los nodos,
#DESP son los desplazamientos que quiero graficar y puede haber términos para evoluciones temporales, por último la 
#dimensionalidad del vector que estoy guardando (desplazamientos son un vector con tres dimensiones)
Desps = gmsh.view.addModelData(desps, 0, 'test', 'NodeData', NodeInfo[0], DESP, numComponents = 3)
Tens = gmsh.view.addModelData(ViewTension, 0, 'test', 'ElementData', ETAGS, Tensiones[:,0].reshape(-1,1), numComponents = 1)

#Si tengo dependencia temporal: (reshape(-1,1) convierte filas en columnas)

#Temps = gmsh.view.add("Temperaturas")
#
#for i in range(100):
#	gmsh.view.addModelData(Temps, i, 'test', 'NodeData', NodeInfo[0], MatrizNodos[:,0].reshape(-1,1)*i, numComponents = 1)

#Guardado de mallado

# gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
# gmsh.write('TestMeshView_in.msh')
# gmsh.view.write(1,"TestMeshView_out.msh", append = True)

#gmsh.model.geo.synchronize()
gmsh.fltk.run()
gmsh.finalize()