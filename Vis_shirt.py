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
MYDEEPBLUE = "#1E90FF"
MYBLACK = "#000000"
MYGRAY = "#C0C0C0"

moisture = True
moi = lambda x:MYBLUE if x else MYGRAY

moisturecolor = moi(moisture)
Z22value = utils.loaddatafile(".calibrate")[1]
slider = 0.5
sliderheight = 0.1
pressfill = [False, True, False]
presscolor = [MYGRAY, MYBLUE, MYGRAY]
pressstyle = ["dashed", "", "dashed"]
presstext = ["", "medium", ""]
#construct figure
plt.rcParams["toolbar"] = "None"
fig = plt.figure(constrained_layout=True, figsize=(12, 6))
fig.canvas.set_window_title("   ")
gs = gridspec.GridSpec(nrows=12, ncols=10, figure=fig)

axZ22 = fig.add_subplot(gs[:, :6])
axSld = fig.add_subplot(gs[:, 6:8])
axMoi = fig.add_subplot(gs[:4, 8:10])
axPre = fig.add_subplot(gs[4:12, 8:10])

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
axZ22.fill_betweenx([1.05, 0], 1e6, 7e6, color=MYBLUE, alpha=0.1)
axZ22.fill_betweenx([1.05, 0], 8e6, 13e6, color=MYBLUE, alpha=0.2)
axZ22.fill_betweenx([1.05, 0], 14e6, 20e6, color=MYBLUE, alpha=0.3)
axZ22.grid(linestyle="-.")

Sld = []
axSld.set_title("slider", y=-0.06)
axSld.spines["left"].set_visible(False)
axSld.spines["right"].set_visible(False)
axSld.spines["top"].set_visible(False)
axSld.spines["bottom"].set_visible(False)
axSld.spines["bottom"].set_position(("data", 0))
axSld.axes.xaxis.set_visible(False)
axSld.axes.yaxis.set_visible(False)
axSld.set_ylim(0, 1.05)
axSld.set_xlim(-0.3, 1.3)
Sld = []
Sld.append(patchs.Rectangle((0.2, 0), 0.6, slider - sliderheight, color=MYGRAY, animated=True))
Sld.append(patchs.Rectangle((0.2, slider - sliderheight), 0.6, 2 * sliderheight, color=MYBLUE, animated=True))
Sld.append(patchs.Rectangle((0.2, slider + sliderheight), 0.6, 1 - slider - sliderheight, color=MYGRAY, animated=True))
for item in Sld:
        axSld.add_patch(item)

serp_height = 0.2
serp_len = 0.5
serp_turn = 4
serp_intv = 0.08
serp_off = 0.02
serpentine_x = []
serpentine_y = []
serpentine_x.append(0.5 - serp_len/2)
serpentine_y.append(serp_height)
for i in range(0, 2 * serp_turn, 2):
        serpentine_x.append(0.5 + serp_len/2)
        serpentine_y.append(serp_height + i * serp_intv)

        serpentine_x.append(0.5 + serp_len/2 + serp_off)
        serpentine_y.append(serp_height + (i + 0.3) * serp_intv)
        serpentine_x.append(0.5 + serp_len/2 + serp_off)
        serpentine_y.append(serp_height + (i + 0.7) * serp_intv)

        serpentine_x.append(0.5 + serp_len/2)
        serpentine_y.append(serp_height + (i + 1) * serp_intv)
        serpentine_x.append(0.5 - serp_len/2)
        serpentine_y.append(serp_height + (i + 1) * serp_intv)

        serpentine_x.append(0.5 - serp_len/2 - serp_off)
        serpentine_y.append(serp_height + (i + 1.3) * serp_intv)
        serpentine_x.append(0.5 - serp_len/2 - serp_off)
        serpentine_y.append(serp_height + (i + 1.7) * serp_intv)

        serpentine_x.append(0.5 - serp_len/2)
        serpentine_y.append(serp_height + (i + 2) * serp_intv)
serpentine_x.append(0.5 + serp_len/2)
serpentine_y.append(serp_height + serp_turn * 2 * serp_intv)

axMoi.set_title("moisture sensor", y=-0.1)
axMoi.spines["left"].set_visible(False)
axMoi.spines["right"].set_visible(False)
axMoi.spines["top"].set_visible(False)
axMoi.spines["bottom"].set_visible(False)
axMoi.axes.xaxis.set_visible(False)
axMoi.axes.yaxis.set_visible(False)
axMoi.set_xlim(0, 1)
axMoi.set_ylim(0, 1)
Moi, = axMoi.plot(serpentine_x, serpentine_y, solid_capstyle="round", linewidth=6, color=moisturecolor, animated=True)

axPre.set_title("pressure sensor", y=-0.09)
axPre.spines["left"].set_visible(False)
axPre.spines["right"].set_visible(False)
axPre.spines["top"].set_visible(False)
axPre.spines["bottom"].set_visible(False)
axPre.axes.xaxis.set_visible(False)
axPre.axes.yaxis.set_visible(False)
axPre.set_xlim(0, 1)
axPre.set_ylim(0, 1.2)
Pre = []
Pre.append(patchs.FancyBboxPatch((0.25, 0.05), 0.5, 0.3, boxstyle="round, pad=0.01", fill=pressfill[0], color=presscolor[0], linewidth=3, linestyle=pressstyle[0], animated=True))
Pre.append(patchs.FancyBboxPatch((0.25, 0.4), 0.5, 0.3, boxstyle="round, pad=0.01", fill=pressfill[1], color=presscolor[1], linewidth=3, linestyle=pressstyle[1], animated=True))
Pre.append(patchs.FancyBboxPatch((0.25, 0.75), 0.5, 0.3, boxstyle="round, pad=0.01", fill=pressfill[2], color=presscolor[2], linewidth=3, linestyle=pressstyle[2], animated=True))
for item in Pre:
        axPre.add_patch(item)
Pre.append(axPre.text(0.18, 0.55, presstext[1], va="center", ha="center", rotation=90, font={"weight":1000, "size":10}))

# construct end

guessvalue = []

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
                print(guessvalue, end="\r")
                cvalue = guessvalue[0]
                moid = False
                slider = 0.1
                sliderheight = 0
                weight = 0
                if cvalue[0] > 100e-12:
                       moid = True
                yield moid, slider, sliderheight, weight, data

def update(data):
        #update value
        global Sld, Moi, Pre, Z22
        global bt, freqlist
        moid, slider, sliderheight, weight, Z22value = data
        colorMoi = MYBLUE if moid == True else MYGRAY
        pressfill = [False, False, False]
        presscolor = [MYGRAY, MYGRAY, MYGRAY]
        pressstyle = ["dashed", "dashed", "dashed"]
        pressheight = [0.2, 0.55, 0.9]
        presstext = ["light", "medium", "heavy"]

        pressfill[weight] = True
        presscolor[weight] = MYBLUE
        pressstyle[weight] = ""        

        Z22.set_data(freqlist, Z22value)
        Sld[0].set_height(slider - sliderheight)     #bar1
        Sld[1].set_y(slider - sliderheight)          #bar2
        Sld[2].set_y(slider + sliderheight)          #bar3
        Sld[2].set_height(1 - slider - sliderheight)
        
        for i in range(3):
                Pre[i].set(fill=pressfill[i], color=presscolor[i], linestyle=pressstyle[i])
        Pre[3].set_position((0.18, pressheight[weight]))
        Pre[3].set_text(presstext[weight])

        Moi.set_color(colorMoi)

        return [*Sld, *Pre, Moi, Z22]

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
T.prepare()
param = {"lxfmr1":C.Lxfmr1, "lxfmr2":2.996e-6, "mask":[[1, 1, 1, 1], [10.6e-6, 4.811e-6, 2.699e-6, 0]]}
T.setparams(param)
input("start")
ani = FuncAnimation(fig, update, data_fetch, interval=150, blit=True)
plt.show()