import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patchs
from matplotlib.animation import FuncAnimation
from triplestagesolver import Tsolver, Csolver
from nanovna import NanoVNA
from time import sleep
import threading
import numpy as np
import math
import random
import utils
import inspect

freqlist = utils.loaddatafile(".calibrate")[0]

MYRED = "#FF6347"
MYBLUE = "#87CEFA"
MYDEEPBLUE = "#1E90FF"
MYBLACK = "#000000"
MYGRAY = "#C0C0C0"

Z22value = utils.loaddatafile(".calibrate")[1]

angle_t = 50
angle_i = 25
angle_m = 25
#construct figure
fig = plt.figure(constrained_layout=True, figsize=(12, 6))
gs = gridspec.GridSpec(nrows=1, ncols=5, figure=fig)

axZ22 = fig.add_subplot(gs[:, :3])
axFig = fig.add_subplot(gs[:, 3:], projection='3d')

axFig.set_title("finger sensors")
axFig.set_xlim(0, 1)
axFig.set_xticks([0, 0.2, 0.4, 0.6, 0.8])
axFig.set_xticklabels(['', 'middle finger', '', 'index finger', ''])
axFig.set_ylim(0, 1)
axFig.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
axFig.set_yticklabels(['', '', '', 'thumb', ''])
axFig.set_zlim(0, 1)
axFig.set_zticks([0, 0.2, 0.4, 0.6, 0.8])
axFig.set_zticklabels([])

mid_r = (angle_m)*math.pi/180
idx_r = (angle_i)*math.pi/180
thb_r = angle_t * math.pi/180
len1 = 0.3
len3 = 0.4

pt_T = [0.9 - len1 * math.sin(thb_r), 0.6 - len1 * math.cos(thb_r), 0.2]
pt_I = [0.7, 0.3 - len3 * math.cos(idx_r), 0.8 - len3 * math.sin(idx_r)]
pt_M = [0.3, 0.3 - len3 * math.cos(mid_r), 0.8 - len3 * math.sin(mid_r)]

pt_00 = [0.7, 1, 0.8]
pt_01 = [0.3, 1, 0.8]
pt_1 = [0.8, 1, 0.4]
pt_t = [0.9, 0.6, 0.2]
pt_i = [0.7, 0.3, 0.8]
pt_m = [0.3, 0.3, 0.8]

ptstr = [pt_i, pt_00, pt_01, pt_m]
ptlst = [pt_00, pt_01, pt_1, pt_t, pt_i, pt_m]
iterptlst = [pt_t, pt_T, pt_i, pt_I, pt_m, pt_M]

Agl = []
axFig.plot3D([pt_00[0], pt_1[0], pt_t[0]], [pt_00[1], pt_1[1], pt_t[1]], [pt_00[2], pt_1[2], pt_t[2]], linewidth=15, solid_capstyle="round", color=MYBLUE)
axFig.plot3D([item[0] for item in ptstr], [item[1] for item in ptstr], [item[2] for item in ptstr], linewidth=15, solid_capstyle="round", color=MYBLUE)
Agl.append(axFig.plot3D([pt_t[0], pt_T[0]], [pt_t[1], pt_T[1]], [pt_t[2], pt_T[2]], linewidth=20, solid_capstyle="round", color=MYRED, animated=True)[0])
Agl.append(axFig.plot3D([pt_i[0], pt_I[0]], [pt_i[1], pt_I[1]], [pt_i[2], pt_I[2]], linewidth=20, solid_capstyle="round", color=MYRED, animated=True)[0])
Agl.append(axFig.plot3D([pt_m[0], pt_M[0]], [pt_m[1], pt_M[1]], [pt_m[2], pt_M[2]], linewidth=20, solid_capstyle="round", color=MYRED, animated=True)[0])
axFig.scatter3D([item[0] for item in ptlst], [item[1] for item in ptlst], [item[2] for item in ptlst], s=100, color=MYBLACK)

axZ22.spines["right"].set_visible(False)
axZ22.spines["top"].set_visible(False)
axZ22.spines["left"].set_position(("data", 0))
axZ22.xaxis.set_ticks_position("bottom")
Z22, = axZ22.plot(freqlist, np.absolute(Z22value))
axZ22.set_xticks([*np.arange(0, 41e6, 5e6), 1e6])
axZ22.set_xticklabels([*np.arange(0, 41, 5), 1])
axZ22.set_ylim(0, 1.05)
axZ22.set_yticks(np.arange(0, 1.05, 0.2))
axZ22.set_title("NanoVNA monitor")
axZ22.set_xlabel("frequency(MHz)")
axZ22.set_ylabel("S11 parameter")
axZ22.fill_betweenx([1.05, 0], 5e6, 8e6, color=MYBLUE, alpha=0.1)
axZ22.fill_betweenx([1.05, 0], 9e6, 12e6, color=MYBLUE, alpha=0.2)
axZ22.fill_betweenx([1.05, 0], 15e6, 20e6, color=MYBLUE, alpha=0.3)
axZ22.grid(linestyle="-.")

def data_fetch():
    global nv, C, T, guessvalue
    while True:
        try:
            nv.fetch_frequencies()
            data = nv.data(0)
            unpacked = C.unpack(data)
            T.loaddata(unpacked)
            T.initiateguess()
            res = T.solve(True)
            guessvalue = res
        except:
            pass
        # print(guessvalue, end="\r")
        cvalue = guessvalue[0]
        angle_t = 10
        angle_i = 10
        angle_m = 10
        if cvalue[0] > 1.55e-11:
             angle_t = 90
        if cvalue[1] > 1.79e-11:
             angle_i = 90
        if cvalue[2] > 2.1e-11:
             angle_m = 90
        # if cvalue[0] > 1.97e-11:
        #      angle_t = 90
        # if cvalue[1] > 1.85e-11:
        #      angle_i = 90
        # if cvalue[2] > 2.58e-11:
        #      angle_m = 90
            
        # if cvalue[0] > 1.65e-11:
        #     angle_t = 90
        # if cvalue[1] > 3.0e-11:
        #     angle_i = 45
        # if cvalue[1] > 3.1e-11:
        #     angle_i = 90
        # if cvalue[1] > 3.2e-11:
        #     angle_i = 150
        # if cvalue[2] > 1.9e-11:
        #     angle_m = 45
        # if cvalue[2] > 2e-11:
        #     angle_m = 90
        # if cvalue[2] > 2.2e-11:
        #     angle_m = 150
        # if cvalue[2] > :
        print(cvalue[0], cvalue[1], cvalue[2])
        yield angle_t, angle_i, angle_m, data

guessvalue = []

def update(data):
        #update value
        global Agl, Z22, pt_t, pt_i, pt_m
        global bt, freqlist, len1, len3
        # angle_m = random.uniform(10, 100)
        # angle_i = random.uniform(10, 100)
        # angle_t = random.uniform(10, 100)
        angle_t, angle_i, angle_m, Z22value = data
        Z22.set_data(freqlist, Z22value)
        
        mid_r = (angle_m)*math.pi/180
        idx_r = (angle_i)*math.pi/180
        thb_r = angle_t * math.pi/180
        pt_T = [0.9 - len1 * math.sin(thb_r), 0.6 - len1 * math.cos(thb_r), 0.2]
        pt_I = [0.7, 0.3 - len3 * math.cos(idx_r), 0.8 - len3 * math.sin(idx_r)]
        pt_M = [0.3, 0.3 - len3 * math.cos(mid_r), 0.8 - len3 * math.sin(mid_r)]
        Agl[0].set_data_3d([pt_t[0], pt_T[0]], [pt_t[1], pt_T[1]], [pt_t[2], pt_T[2]])
        Agl[1].set_data_3d([pt_i[0], pt_I[0]], [pt_i[1], pt_I[1]], [pt_i[2], pt_I[2]])
        Agl[2].set_data_3d([pt_m[0], pt_M[0]], [pt_m[1], pt_M[1]], [pt_m[2], pt_M[2]])

        return [*Agl, Z22]

print("initialize")
nv = NanoVNA()
nv.set_sweep(1e6, 40e6)
nv.set_frequencies(1e6, 40e6)

C = Csolver()
T = Tsolver()
input("calibrate")
nv.fetch_frequencies()
data = nv.data(0)
C.calibrate(data)
# print(C.Csma, C.Lsma, C.Lxfmr1, C.Rl)
print()
param = {"lxfmr1":C.Lxfmr1, "lxfmr2":2.879e-6, "mask":[[1, 1, 1, 1], [22e-6, 10e-6, 1e-6, 0]]}
T.prepare()
T.setparams(param)
input("start")
ani = FuncAnimation(fig, update, data_fetch, interval=150, blit=True)
plt.show()