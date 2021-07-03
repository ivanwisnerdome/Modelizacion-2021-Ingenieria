import numpy as np
import pdb

def teo(x):
    """Solución teórica de la ecuación de la guía"""
    return (4/1.3)*( np.exp(0.8*x)-np.exp(-0.5*x) ) + 2*np.exp(-0.5*x)

def F(_x, _y):
    """Ecuación Diferencial"""
    return 4*np.exp(0.8*_x) - 0.5 * _y

def paso_euler(_dx, _xo, _yo, _F):
    """
    Solución de un paso por el método de Euler
    Parámetros:
    ===========
    _dx = tamaño de paso
    _xo = posición inicial
    _yo = condición inicial
    _F = función característica de la ecuación diferencial
    
    Retorna
    ===========
    valor de la solución en x+dx
    """
    return _yo + _dx*_F(_xo, _yo)

def paso_RK(_dx, _xo, _yo, _F):
    """
    Solución de un paso por el método de Runge - Kutta de orden 4
    Parámetos:
    ===========
    _dx = tamaño de paso
    _xo = posición inicial
    _yo = condición inicial
    _F = función característica de la ecuación diferencial
    
    Retorna
    ==========
    valor de la solución en x+dx
    """
    K1 = _F(_xo, _yo)
    K2 = _F(_xo + 0.5*_dx, _yo+0.5*K1*_dx)
    K3 = _F(_xo + 0.5*_dx, _yo+0.5*K2*_dx)
    K4 = _F(_xo +     _dx, _yo+    K3*_dx)
    return _yo + (1/6)*(K1+2*K2+2*K3+K4)*_dx

def dfsolve(dx=0.1, xo = 0, xf = 4, yo = 2, func=F, method=paso_euler):
    """
    Motor de solución para una ecuación temporal para todo el intervalo del problema
    
    Parámetros:
    ==========
    dx: tamaño de paso, = 0.1
    xo: posición inicial, =0
    xf: posición final (límite del intervalo), = 4
    yo: condición inicial, =2
    func: ecuación característica de la ecuacion diferencial (callable) = F
    method: método de solución(callable) = paso_euler
    
    Retorna:
    ==========
    X: vector de las posiciones donde se resolvió la ecuación
    SOL: valores de Y con la solución encontrada
    """
    X = np.linspace(xo, xf, int((xf-xo)/dx+1))
    SOL = [yo]
    for i, x in enumerate(X[:-1]):
        SOL.append(
            method(dx, x, SOL[-1], func )
        )
    return X, SOL
        