import numpy as np
import matplotlib.pyplot as plt

S0 = 0.9
I0 = 0.1
R0 = 0.0

def funcion (i):
       #MDF-COMMENT lo que no me gusta de este approach es que necesitas 3 llamadas a la función ,
       #MDF-COMMENT cuando podrías hacerlo con una sola
       #MDF-COMMENT y por otro lado, los valores del paso anterior (S,I,R) tendrían que ser los verdaderos 
       #MDF-COMMENT argumentos de tu función.
	dS = lambda t,S,I,R: -0.35*I*S
	dI = lambda t,S,I,R: 0.35*I*S - 0.1*I
	dR = lambda t,S,I,R: 0.1*I
	Der = [dS, dI, dR]

        #MDF-COMMENT lo que estar retornando es la función !
	return Der[i]

#MDF-COMMENT a ver si me sale.
def funcionv(Y, _gamma, _beta):
#MDF-COMMENT a demas te conviene poner como argumentos gamma y beta, que son los parametros del modelo
    """
    Parámetros
    ==========
    Y = [S, I, R]
    """
    dY = np.array( [ -_gamma*Y[1]*Y[0],  _gamma*Y[1]*Y[0]-_beta*Y[1],  _beta*Y[1]  ]).reshape(-1,1)
    return dY



t = np.linspace(0, 100, 1000)
k = np.zeros(4)
Y = np.zeros((3,len(t)))
Y[0,0] = S0
Y[1,0] = I0
Y[2,0] = R0

Espaciado = t[1] - t[0]
#MDF-COMMENT 
Espaciado = np.diff(t)

gamma = 0.5
beta = 0.01
for i in range(1,len(t)-1):
#    for j in range(0,3):
    K1 = funcionv(Y[:, i-1], gamma, beta)
    K2 = funcionv(Y[:, i-1].reshape(-1,1)+K1*Espaciado[i]*0.5, gamma, beta)
    K3 = funcionv( Y[:, i-1].reshape(-1,1) +K2*Espaciado[i]*0.5, gamma, beta)
    K4 = funcionv( Y[:, i-1].reshape(-1,1) +K3*Espaciado[i],  gamma, beta)
    Y[:, i] = (Y[:, i-1].reshape(-1,1) + (1/6)*Espaciado[i]*(K1+2*K2*2*K3+K4)).ravel()
#MDF-COMMENT lo que habías hecho no estaba mal, pero esto de arriva es un poco menos dificil de leer,
#MDF-COMMENT lo unico que hace un poco de lío es el tema de que lo definiste todo como vector columna
#            G = funcionv(Y[:,-1])
#MDF-COMMENT            G = funcion(j)
#MDF-COMMENT            k[0] = G(t[i-1], Y[0,i-1], Y[1,i-1], Y[2,i-1])
#MDF-COMMENT            k[1] = G(t[i-1] + Espaciado/2 , Y[0,i-1] + k[0]*Espaciado/2, Y[1,i-1] + k[0]*Espaciado/2, Y[2,i-1] + k[0]*Espaciado/2)
#MDF-COMMENT            k[2] = G(t[i-1] + Espaciado/2 , Y[0,i-1] + k[1]*Espaciado/2, Y[1,i-1] + k[1]*Espaciado/2, Y[2,i-1] + k[1]*Espaciado/2)
#MDF-COMMENT            k[3] = G(t[i-1] + Espaciado, Y[0,i-1] +k[2]*Espaciado, Y[1,i-1] +k[2]*Espaciado, Y[2,i-1] +k[2]*Espaciado)
#MDF-COMMENT            Y[j,i] = Y[j,i-1] + 1/6*(k[0] + 2*k[1] + 2*k[2] + k[3])*Espaciado

plt.plot(t,Y[0,:])
plt.plot(t,Y[1,:])
plt.plot(t,Y[2,:])
#MDF-COMMENT  plt.show()
#MDF-COMMENT a demas de mostrar el grafico, podes guardarlo en un archivo
