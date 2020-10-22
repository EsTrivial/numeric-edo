from numpy import *
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# impresión de resultados
"""
printSoln(X,Y,freq,[exact]).
Imprime los valores de X e Y cada 'freq' pasos
freq = n imprime cada n pasos
freq = 0 imprime solo el paso final
Imprime (opcionalmente) el error entre la solución 
exacta (exact) y la calculada (Y)
Retorna (opcionalmente) el máximo valor absoluto del error
"""
def printSoln(X, Y, freq, exact = []):
    # imprime cabecera
    def printHead(n):
        print("\n", end="")
        print("{:^15}".format('x'), end="")
        if n == 1:
            print("{:^15}".format("y"), end="")
        else:
            for i in range (n):
                print("{:^15}".format("y[{}]".format(i)), end="")
        if len(exact):
            print("{:^14}".format("Error"),end="")
        print()
    # impresión línea a línea
    def printLine(x,y,n,z = None):
        print("{:14.10f}".format(x),end="")
        for i in range (n):
            print("{:14.10f}".format(y[i]),end="")
        if z != None:
            print("{:14.10f}".format(z),end="")
        print()
    m = len(Y)
    try: 
        n = len(Y[0])
    except TypeError: 
        n = 1
    if freq == 0: 
        freq = m
    printHead(n)
    if len(exact):
        exact = abs(exact-Y[:,0])
        ret = exact.max()
    else:
        exact = m*[None]
        ret = None
    for i in range(0,m,freq):
            printLine(X[i],Y[i],n,exact[i])
    if i != m - 1: 
        printLine(X[m - 1],Y[m - 1],n,exact[m-1])
    return ret
        
# método de Euler
""" 
X,Y = euler(F,x,y,xStop,h).
Método de Euler para resolver un problema 
de valores iniciales y' = F(x,y), donde
y = (y[0],y[1],...y[n-1]).
    x,y = condiciones iniciales
    xStop = valor final para x
    h = paso de integración
    F = función escrita por el usuario que retorna un 
    array F(y,x) = (y'[0],y'[1],...,y'[n-1]).
"""
def euler(F, x, y, xStop = 1., h = 0.1):
    if x >= xStop:
        print("Error: no hay intervalo. Calculando en (0,1)")
        x = 0.
        xStop = 1.
    X = arange(x,xStop+h,h)
    Y = zeros((X.shape[0],y.shape[0]))
    Y[0,:] = y
    for i in range(X.shape[0]-1):
        Y[i+1,:] = Y[i,:] + h*F(Y[i,:],X[i])
    return X,Y


# método de Euler Modificado
""" 
X,Y = eulerModificado(F,x,y,xStop,h).
Método de Euler Modificado para resolver un problema 
de valores iniciales y' = F(x,y), donde
y = (y[0],y[1],...y[n-1]).
    x,y = condiciones iniciales
    xStop = valor final para x
    h = paso de integración
    F = función escrita por el usuario que retorna un 
    array F(y,x) = (y'[0],y'[1],...,y'[n-1]).
"""
def eulerModificado(F, x, y, xStop = 1., h = 0.1):
    if x >= xStop:
        print("Error: no hay intervalo. Calculando en (0,1)")
        x = 0.
        xStop = 1.
    X = arange(x,xStop+h,h)
    Y = zeros((X.shape[0],y.shape[0]))
    Y[0,:] = y
    for i in range(X.shape[0]-1):
        Y[i+1,:] = Y[i,:] + h/2*( F(Y[i,:],X[i]) + F(Y[i,:]+h*F(Y[i,:],X[i]),X[i+1]))
    return X,Y


# método de Medio Paso
""" 
X,Y = medioPaso(F,x,y,xStop,h).
Método de Euler Modificado para resolver un problema 
de valores iniciales y' = F(x,y), donde
y = (y[0],y[1],...y[n-1]).
    x,y = condiciones iniciales
    xStop = valor final para x
    h = paso de integración
    F = función escrita por el usuario que retorna un 
    array F(y,x) = (y'[0],y'[1],...,y'[n-1]).
"""
def medioPaso(F, x, y, xStop = 1., h = 0.1):
    if x >= xStop:
        print("Error: no hay intervalo. Calculando en (0,1)")
        x = 0.
        xStop = 1.
    X = arange(x,xStop+h,h)
    Y = zeros((X.shape[0],y.shape[0]))
    Y[0,:] = y
    for i in range(X.shape[0]-1):
        Y[i+1,:] = Y[i,:] + h*(F(Y[i,:] + h/2*F(Y[i,:],X[i]),X[i]+h/2))
    return X,Y



# integración mediante odeint de scipy
""" 
X,Y = odesolve(F,x,y,xStop,h).
Método de Euler Modificado para resolver un problema 
de valores iniciales y' = F(x,y), donde
y = (y[0],y[1],...y[n-1]).
    x,y = condiciones iniciales
    xStop = valor final para x
    h = paso de integración
    F = función escrita por el usuario que retorna un 
    array F(y,x) = (y'[0],y'[1],...,y'[n-1]).
"""
def odesolve(F, x, y, xStop = 1., h = 0.1):
    X = arange(x,xStop+h,h)
    Y = odeint(F, y, X)
    return X,Y



# dibujo de resultados
"""
dibujar(X,Y,[yExact])
Dibuja la solución la primera componente de la 
solución del sistema, o la solución de la edo
Dibuja opcionalmente la solución exacta
"""
def dibujar(X,*Y):
    if len(Y)>1:
        for i in range(len(Y)-1):
            plt.plot(X,Y[i][:,0],'o',label="Num "+str(i+1))
        plt.plot(X,Y[len(Y)-1],'-',label="Exacta")
        plt.legend(('Numerica','Exacta'),loc=0)
    else:
        plt.plot(X,Y[0][:,0],'-o',label = "Numérica")
    plt.legend(loc=0)
    plt.xlabel('x'); plt.ylabel('y')
