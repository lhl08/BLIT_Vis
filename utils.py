from sympy import *
import numpy as np
import matplotlib.pyplot as plt

import warnings

# paramdict = {
#     "purse": 
#         {"lxfmr1":C.Lxfmr1, "lxfmr2":3.139e-6, "mask":[[1, 1, 0, 1], [6.12e-6, 6.62e-6, 15e-12, 0]]},
#     "garment":
#         {"lxfmr1":C.Lxfmr1, "lxfmr2":2.996e-6, "mask":[[1, 1, 1, 1], [10.6e-6, 4.811e-6, 2.699e-6, 0]]},
#     "phonecase":
#         {"lxfmr1":C.Lxfmr1, "lxfmr2":2.886e-6, "mask":[[1, 1, 0, 1], [8.87e-6, 5.92e-6, 13.05e-12, 0]]},
#     "glove":
#         {"lxfmr1":C.Lxfmr1, "lxfmr2":2.879e-6, "mask":[[1, 1, 1, 1], [12.594e-6, 6.041e-6, 2.726e-6, 0]]}
# }

warnings.filterwarnings("ignore")

Z_o = symbols('Z_o', real=True)
S = symbols('S')
Zi = symbols('Zi')

Z = (Z_o * (1 + S))/(1 - S)
Z = Z.subs(Z_o, 50)
Zfunc = lambdify((S), Z, "numpy")

S11 = (Zi - Z_o)/(Zi + Z_o)
S11 = S11.subs(Z_o, 50)
Sfunc = lambdify((Zi), S11, "numpy")

def loaddatafile(name):
    arx = []
    arreal = []
    arimag = []
    with open(name + ".data", mode="r") as datafile:
        line = datafile.readline()
        while line:
            lst = line.strip().split(" ")
            lst = [float(lst[0]), complex(lst[1])]
            arx.append(lst[0])
            arreal.append(lst[1].real)
            arimag.append(lst[1].imag)

            line = datafile.readline()

        return np.array(arx), np.array(arreal) + 1j * np.array(arimag)

def loadfile(name):
    arx = []
    arreal = []
    arimag = []
    with open(name + ".csv", mode="r") as csvfile:
        line = csvfile.readline()
        while line:
            lst = line.strip().split(";")
            lst = [item.strip() for item in lst]
            if float(lst[0]) > 40e6:
                break
            arx.append(float(lst[0]))
            arreal.append(float(lst[1]))
            arimag.append(float(lst[2]))
            line = csvfile.readline()
    return np.array(arx), np.array(arreal) + 1j * np.array(arimag)

def diff_2(y, x):
    return ((y[2] - y[0])/(x[2] - x[0]) - (y[1] - y[0])/(x[1] - x[0]))/(x[2] - x[1])

def threeptsfit(y, x):
    #to get x0, y0, and a -> y = a(x - x0)**2 + y0
    k = diff_2(y, x)
    a = (y[1] - y[0])/(x[1] - x[0])
    xm = (x[0] + x[1] - a/k)/2
    ym = y[0] + (xm - x[0])*(k*(xm - x[1]) + a)
    return k, xm, ym

def drawcomplexplot(arx, arZlist, filename="test", ylabel="test"):
    ax = plt.gca()
    ax.spines["top"].set_color("None")
    ax.spines["right"].set_color("None")
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    ax.spines["left"].set_position(("data", 0))
    ax.grid(linestyle="dashed")
    ax.set_ylabel(ylabel)
    for i in range(len(arZlist)):
        plt.plot(arx, np.real(arZlist[i]), label="real-"+str(i), linewidth="1")
        plt.plot(arx, np.imag(arZlist[i]), label="imag-"+str(i), linewidth="1")
    plt.legend()
    plt.savefig(filename + ".jpg", dpi=100)
    plt.clf()

def drawplot(arx, arylist, filename="test", ylabel="test"):
    ax = plt.gca()
    ax.spines["top"].set_color("None")
    ax.spines["right"].set_color("None")
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    ax.spines["left"].set_position(("data", 0))
    ax.grid(linestyle="dashed")
    ax.set_xlabel("frequency(MHz)")
    ax.set_ylabel(ylabel)
    for ary in arylist:
        plt.plot(arx, ary, label=ylabel, linewidth="1")
    plt.legend()
    plt.savefig(filename + ".jpg", dpi=800)
    plt.clf()

