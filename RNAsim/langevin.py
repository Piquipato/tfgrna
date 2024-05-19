import matplotlib.pyplot as plt
import numpy as np

def lange(x):
    return (1/np.tanh(x)) - (1/x)

def inv_lange(x):
    return (x*(3-1.00651*x**2-0.962251*x**4+1.47353*x**6-0.46953*x**8))/((1-x)*(1+1.01524*x))

def lnZ(x):
    return np.log(4*np.pi*np.sinh(x)) - np.log(x)

if __name__ == '__main__':
    fig, ax = plt.subplots()
    x = np.linspace(0, 1, 500)
    ax.plot(x, -300*lnZ(inv_lange(x)))
    plt.show()