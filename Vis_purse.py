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

freqlist = utils.loaddatafile(".calibrate")[0]

MYRED = "#FF6347"
MYBLUE = "#87CEFA"
MYBLACK = "#000000"
MYGRAY = "#C0C0C0"

button = False
bt = lambda x:MYBLUE if x else MYGRAY

obj = "coins"
# obj = "_____"
colorbtn = bt(button)
stretch = 0.2
stretchcolor = MYBLUE
Z22value = utils.loaddatafile(".calibrate")[1]

#construct figure
plt.rcParams["toolbar"] = "None"
fig = plt.figure(constrained_layout=True, figsize=(12, 6))
fig.canvas.set_window_title("   ")
gs = gridspec.GridSpec(ncols=4, nrows=3, figure=fig)

axZ22 = fig.add_subplot(gs[:, :3])
axObj = fig.add_subplot(gs[2, 3])
axBtn = fig.add_subplot(gs[1, 3])
axSch = fig.add_subplot(gs[0, 3])

axZ22.spines["right"].set_visible(False)
axZ22.spines["top"].set_visible(False)
axZ22.spines["left"].set_position(("data", 0))
axZ22.xaxis.set_ticks_position("bottom")
Z22, = axZ22.plot(freqlist, np.absolute(Z22value), animated=True)
axZ22.set_xticks([*np.arange(0, 41e6, 5e6), 1e6])
axZ22.set_xticklabels([*np.arange(0, 41, 5), 1])
axZ22.set_ylim(0, 1.05)
axZ22.set_yticks(np.arange(0, 1.05, 0.2))
axZ22.set_title("NanoVNA monitor")
axZ22.set_xlabel("frequency(MHz)")
axZ22.set_ylabel("S11 parameter")
axZ22.fill_betweenx([1.05, 0], 6e6, 8e6, color=MYBLUE, alpha=0.1)
axZ22.fill_betweenx([1.05, 0], 9e6, 13e6, color=MYBLUE, alpha=0.2)
axZ22.fill_betweenx([1.05, 0], 15e6, 23e6, color=MYBLUE, alpha=0.3)
axZ22.grid(linestyle="-.")

axBtn.set_title("button", y=-0.24)
axBtn.spines["left"].set_visible(False)
axBtn.spines["right"].set_visible(False)
axBtn.spines["top"].set_visible(False)
axBtn.spines["bottom"].set_visible(False)
axBtn.axes.xaxis.set_visible(False)
axBtn.axes.yaxis.set_visible(False)
axBtn.set_aspect("equal")
cirbtn = plt.Circle((0.5, 0.5), 0.4, fill=True, color=colorbtn, animated=True)
axBtn.add_artist(cirbtn)
Btn = cirbtn

axObj.spines["left"].set_visible(False)
axObj.spines["right"].set_visible(False)
axObj.spines["top"].set_visible(False)
axObj.spines["bottom"].set_visible(False)
axObj.axes.xaxis.set_visible(False)
axObj.axes.yaxis.set_visible(False)
axObj.set_xlim(0, 10)
axObj.set_ylim(0, 10)
axObj.set_title("object detector", y=-0.1)
axObj.add_patch(patchs.FancyBboxPatch((1, 2), 8, 6, boxstyle="round, pad=0.5", fill=False, linewidth=2, edgecolor=MYBLACK))
textobj = axObj.text(5, 5, obj, size=40, style="normal", weight=500, ha="center", va="center", animated=True)
Obj = textobj

Sch = []
axSch_height = 0.4
axSch_floater = 0.03
axSch.spines["left"].set_visible(False)
axSch.spines["right"].set_visible(False)
axSch.spines["top"].set_visible(False)
axSch.spines["bottom"].set_visible(False)
axSch.axes.xaxis.set_visible(False)
axSch.axes.yaxis.set_visible(False)
axSch.set_xlim(-0.1, 1.1)
axSch.set_ylim(0, 1)
axSch.set_title("stretch sensor", y=-0.2)
axSch.plot([0, 1], [axSch_height, axSch_height], linewidth=5, linestyle="dashdot", color=MYBLACK)
axSch.plot([1, 1], [axSch_height + axSch_floater, axSch_height - axSch_floater], linewidth=5, color=MYBLACK)
Sch.append(axSch.plot([0, 0], [axSch_height + axSch_floater, axSch_height - axSch_floater], linewidth=10, color=stretchcolor, solid_capstyle="round", animated=True)[0])
Sch.append(axSch.plot([0, stretch], [axSch_height, axSch_height], linewidth=8, color=stretchcolor, solid_capstyle="round", animated=True)[0])
Sch.append(axSch.plot([stretch, stretch], [axSch_height + axSch_floater, axSch_height - axSch_floater], linewidth=10, color=stretchcolor, solid_capstyle="round", animated=True)[0])

guessvalue = []

def data_fetch():
    global nv, C, T, guessvalue
    while True:
        try:
            nv.fetch_frequencies()
            data = nv.data(0)
            unpacked = C.unpack(data)
            T.loaddata(unpacked)
            T.initiateguess(jump0=True)
            res = T.solve(True, jump0=True)
            guessvalue = res
        except:
            pass
        print(guessvalue, end="\r")

        #add logic
        colorbtn = MYGRAY
        colorsch = MYBLUE
        stretch = 0.3
        # text = random.choice(["coin", "keys"])
        text = ""
        if guessvalue[0][2] > 5.4e-6:
            text = "ID card"
        if guessvalue[0][2] < 4.8e-6:
            text = "key"
        if guessvalue[0][2] < 4.2e-6:
            text = "coins"
        
        # if guessvalue[0][1] > 3e-11:
        #     colorbtn = MYBLUE
        # else:
        #     if guessvalue[1][0] > 0.018:
        #         stretch = 0.95
        #         colorsch = MYRED

        #logic end
        yield np.absolute(data), colorbtn, colorsch, stretch, text

def update(data):
    #update value
    global Sch, Btn, Obj, Z22
    global bt, freqlist
    Z22value, colorbtn, colorsch, stretch, text = data

    Z22.set_data(freqlist, Z22value)
    Btn.set_color(colorbtn)
    Obj.set_text(text)
    Sch[0].set_color(colorsch)
    Sch[1].set_color(colorsch)
    Sch[2].set_color(colorsch)
    Sch[1].set_data([0, stretch], [axSch_height, axSch_height])
    Sch[2].set_data([stretch, stretch], [axSch_height + axSch_floater, axSch_height - axSch_floater])
    return [*Sch, Btn, Obj, Z22]

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
print(C.Csma, C.Lsma, C.Lxfmr1, C.Rl)
param = {"lxfmr1":C.Lxfmr1, "lxfmr2":3.139e-6, "mask":[[1, 1, 0, 1], [6.12e-6, 6.62e-6, 15e-12, 0]]}
T.prepare()
T.setparams(param)
input("start")
ani = FuncAnimation(fig, update, data_fetch, interval=150, blit=True)
plt.show()