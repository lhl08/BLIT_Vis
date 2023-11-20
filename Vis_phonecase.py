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

freqlist = np.linspace(1e6, 40e6, 101)

MYRED = "#FF6347"
MYBLUE = "#87CEFA"
MYBLACK = "#000000"
MYGRAY = "#C0C0C0"

button1 = True
button2 = False
bt = lambda x:MYBLUE if x else MYGRAY
colorBup = bt(button1)
colorBdn = bt(button2)
slider = 0.5
sliderheight = 0.1
angle = 30 * math.pi / 180
Z22value = utils.loaddatafile(".calibrate")[1]

#construct figure
plt.rcParams["toolbar"] = "None"
fig = plt.figure(constrained_layout=True, figsize=(12, 8))
fig.canvas.set_window_title("   ")
gs = gridspec.GridSpec(nrows=12, ncols=12, figure=fig)

axZ22 = fig.add_subplot(gs[1:11, :6])
axSld = fig.add_subplot(gs[:, 6:8])
axBup = fig.add_subplot(gs[:3, 8:10])
axBdn = fig.add_subplot(gs[:3, 10:12])
axAgl = fig.add_subplot(gs[3:12, 8:12])

axSld.set_title("slider", y=-0.06)
axSld.spines["left"].set_visible(False)
axSld.spines["right"].set_visible(False)
axSld.spines["top"].set_visible(False)
axSld.spines["bottom"].set_visible(False)
axSld.spines["bottom"].set_position(("data", 0))
axSld.axes.xaxis.set_visible(False)
axSld.axes.yaxis.set_visible(False)
axSld.set_ylim(0, 1)
axSld.set_xlim(0, 1)
Sld = []
Sld.append(patchs.Rectangle((0.2, 0), 0.6, slider - sliderheight, color=MYGRAY, animated=True))
Sld.append(patchs.Rectangle((0.2, slider - sliderheight), 0.6, 2 * sliderheight, color=MYBLUE, animated=True))
Sld.append(patchs.Rectangle((0.2, slider + sliderheight), 0.6, 1 - slider - sliderheight, color=MYGRAY, animated=True))
for item in Sld:
        axSld.add_patch(item)

axBup.set_title("button (up)", y=-0.24)
axBup.spines["left"].set_visible(False)
axBup.spines["right"].set_visible(False)
axBup.spines["top"].set_visible(False)
axBup.spines["bottom"].set_visible(False)
axBup.axes.xaxis.set_visible(False)
axBup.axes.yaxis.set_visible(False)
axBup.set_aspect("equal")
cirup = plt.Circle((0.5, 0.5), 0.4, fill=True, color=colorBup, animated=True)
axBup.add_artist(cirup)

axBdn.set_title("button (down)", y=-0.24)
axBdn.spines["left"].set_visible(False)
axBdn.spines["right"].set_visible(False)
axBdn.spines["top"].set_visible(False)
axBdn.spines["bottom"].set_visible(False)
axBdn.axes.xaxis.set_visible(False)
axBdn.axes.yaxis.set_visible(False)
axBdn.set_aspect("equal")
cirdn = plt.Circle((0.5, 0.5), 0.4, fill=True, color=colorBdn, animated=True)
axBdn.add_artist(cirdn)
Btn = [cirup, cirdn]

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
axZ22.fill_betweenx([1.05, 0], 5e6, 8e6, color=MYBLUE, alpha=0.1)
axZ22.fill_betweenx([1.05, 0], 10e6, 15e6, color=MYBLUE, alpha=0.2)
axZ22.fill_betweenx([1.05, 0], 19e6, 25e6, color=MYBLUE, alpha=0.3)
axZ22.grid(linestyle="-.")

axAgl.set_title("anglular sensor", y=-0.09)
axAgl.spines["left"].set_visible(False)
axAgl.spines["right"].set_visible(False)
axAgl.spines["top"].set_visible(False)
axAgl.spines["bottom"].set_visible(False)
axAgl.axes.xaxis.set_visible(False)
axAgl.axes.yaxis.set_visible(False)
axAgl.set_ylim(0, 1.5)
axAgl.set_xlim(0, 1.5)
axAgl.set_aspect("equal")
axAgl.add_artist(plt.Circle((0, 0), 1, fill=False, edgecolor="#000000"))
axAgl.plot([0, 1.45], [0, 0], color="#000000", linewidth=5)
axAgl.plot([0, 0], [0, 1.45], color="#000000", linewidth=5, linestyle="dashed")
anglept = (1.45*math.cos(angle), 1.45*math.sin(angle))
Agl, = axAgl.plot([0, anglept[0]], [0, anglept[1]], color="#000000", animated=True) #plot1
ccle = np.linspace(angle, 0, 20)
Aglfill = axAgl.fill_between(np.append(0, np.cos(ccle)), np.append(0, np.sin(ccle)), \
        interpolate=True, color=MYBLUE, alpha=0.8, animated=True) #color
# construct end

def update(data):
        #update value
        global Sld, Agl, Aglfill, axAgl, Btn, Z22
        global bt, freqlist
        Z22value, On_bup, On_bdn, angvalue, slider, sliderheight = data
        colorBup = bt(On_bup)
        colorBdn = bt(On_bdn)
        angle = angvalue * math.pi / 180        

        Z22.set_data(freqlist, Z22value)
        Sld[0].set_height(slider - sliderheight)     #bar1
        Sld[1].set_y(slider - sliderheight)          #bar2
        Sld[2].set_y(slider + sliderheight)          #bar3
        Sld[2].set_height(1 - slider - sliderheight)
        Btn[0].set_color(colorBup)
        Btn[1].set_color(colorBdn)
        anglept = (1.45*math.cos(angle), 1.45*math.sin(angle))
        ccle = np.linspace(angle, 0, 10)
        pth = [[0, 0]]
        for ag in ccle:
                pth.append([math.cos(ag), math.sin(ag)])
        Agl.set_data([0, anglept[0]], [0, anglept[1]])
        Aglfill.set_paths([pth])
        return [*Sld, *Btn, Agl, Aglfill, Z22]

guessvalue = []

def data_fetch():
        global nv, C, T, guessvalue
        while True:
                nv.fetch_frequencies()
                data = nv.data(0)
                unpacked = C.unpack(data)
                try:
                        T.loaddata(unpacked)
                        T.initiateguess()
                        res = T.solve()
                        guessvalue = res
                except:
                        pass
                print(guessvalue, T.true_pkstarget, end="\r")
                # btns
                guess_btn_up = False
                guess_btn_dn = False
                guess_angle = 0
                guess_slider = 0.5
                guess_height = 0.1
                if guessvalue[2] < 3.25e-6:
                        guess_height = 0
                        if guessvalue[2] < 2.9e-6:
                                guess_angle = 0
                        else:
                                guess_angle = (guessvalue[2]*1e6 - 2.9)*200
                else:
                        if guessvalue[1] > 1.55e-11:
                                guess_height = 0
                                if guessvalue[1] > 2.3e-11:
                                        guess_btn_up = True
                                        guess_btn_dn = True
                                elif guessvalue[1] > 1.95e-11:
                                        guess_btn_up = True
                                else:
                                        guess_btn_dn = True
                        else:
                                if guessvalue[0] > 4.4e-11:
                                        guess_slider = 0.9
                                        guess_height = 0.1
                                elif guessvalue[0] > 4.2e-11:
                                        guess_slider = 0.7
                                        guess_height = 0.1
                                elif guessvalue[0] > 4e-11:
                                        guess_slider = 0.5
                                        guess_height = 0.1
                                elif guessvalue[0] > 3.97e-11:
                                        guess_slider = 0.3
                                        guess_height = 0.1
                                elif guessvalue[0] > 3.94e-11:
                                        guess_slider = 0.1
                                        guess_height = 0.1
                                else:
                                        guess_height = 0
                
                yield np.absolute(data), guess_btn_up, guess_btn_dn, \
                        guess_angle, guess_slider, guess_height

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
param = {"lxfmr1":C.Lxfmr1, "lxfmr2":2.886e-6, "mask":[[1, 1, 0, 1], [8.87e-6, 5.92e-6, 13.05e-12, 0]]}
T.prepare()
T.setparams(param)
input("start")
ani = FuncAnimation(fig, update, data_fetch, interval=150, blit=True)
plt.show()