from utils import threeptsfit, Zfunc, diff_2
from nanovna import NanoVNA
from math import pi, sqrt
from scipy.signal import find_peaks
from scipy.optimize import leastsq, least_squares
import numpy as np
import warnings

warnings.filterwarnings("ignore")

class Tsolver(): # use triple-direction iteration to get rlc pairs values
    def __init__(self):
        '''
        clineloss & rmloss : for LM algorithm
        '''
        self.ptsnum = 101
        self.loadnum = 3
        self.freqlist = np.linspace(1e6, 40e6, self.ptsnum)
        self.lxfmr1 = 2.3e-6
        self.lxfmr2 = 2e-6
        self.thresh = int(self.ptsnum * (24 / 40))

        self.true_ttest = None
        self.true_pkstarget = None
        self.true_pksrange = None

        self.guess_mask = [1, 0, 0, 1]
        self.guess_mxfmr = 1e-6
        self.guess_rlist = [0.1, 0.1, 0.1, 0.1]
        self.guess_rmloss = None
        self.guess_rmloss2 = None
        self.guess_rmloss3 = None
        self.guess_Clineloss = None

        self.known_lclist = None
        self.guess_lclist = None

        self.subs_guess_lclist = None
        self.subs_guess_rlist = None    
        self.subs_guess_mxfmr = None    

    def Z22(self, mxfmr, rlist, llist, clist):
        '''
        Z_all = Z_11 + Z_12^2 / Z_22
        Z_all : 总阻抗
        Z_11  : 原电路阻抗
        Z_12  : 互感部分阻抗
        Z_22  : 副电路阻抗
        t = Z_12^2 / Z_22
        '''
        yrlc = 0
        for i in range(self.loadnum + 1):
            yrlc += 1/(rlist[i] + 2j * pi * self.freqlist * llist[i] + 1/(2j * pi * self.freqlist * clist[i]))
        t = (4 * pi * pi * self.freqlist**2 * mxfmr**2) / (1/yrlc + 2j * pi * self.freqlist * self.lxfmr2)
        return np.array(t)

    def loaddata(self, t):
        self.true_ttest = t
        pks = list(find_peaks(t, distance=3)[0])
        maxpks = pks[:self.loadnum]
        pksheight = [t[i].real for i in maxpks]
        for item in pks[self.loadnum:]:
            if item > self.thresh:
                break
            height = t[item].real
            minheight = min(pksheight)
            minindex = pksheight.index(minheight)
            if minheight < height:
                del pksheight[minindex]
                del maxpks[minindex]
                maxpks.append(item)
                pksheight.append(height)
        maxpks = pks[:self.loadnum]
        maxpks = np.array(maxpks)
        pksfit = []
        for pk in maxpks[:self.loadnum]:
            x = [pk - 1, pk, pk + 1]
            y = np.real([t[x[0]], t[x[1]], t[x[2]]])
            pksfit.append(threeptsfit(y, x)[1])
        cpk = self.ptsnum - 1
        ret = 1
        if pks[-1] > self.thresh and t[pks[-1]].real > t[-1].real:
            ret = 2
            maxZ = 0
            maxfreq = -1
            for freq in pks[3:]:
                if t[freq] > maxZ:
                    maxfreq = freq
                    maxZ = t[freq]
            cpk = maxfreq
            cpk = pks[-1] ####
            x = [cpk - 1, cpk, cpk + 1]
            y = np.real([t[x[0]], t[x[1]], t[x[2]]])
            cpk = threeptsfit(y, x)[1]
        pksfit.append(cpk)
        pkstot = np.array(pksfit)
        self.true_pkstarget = pkstot
        # print(self.true_pkstarget)
        return ret

    def initiateguess(self, jump0=False):
        self.guess_mxfmr = 0.2 * sqrt(self.lxfmr2 * self.lxfmr1)
        self.guess_rlist = [20, 20, 5, 10]
        self.guess_lclist = []
        for i in range(self.loadnum + 1):
            if jump0 and i == 0:
                self.guess_lclist.append(54.45e-12)
            else:
                if self.guess_mask[i] == 1: #know l, guess c
                    self.guess_lclist.append(max(1/((2*pi*(1e6 + (40e6 - 1e6) * self.true_pkstarget[i]/(self.ptsnum - 1)))**2 * (self.known_lclist[i] + self.lxfmr2)), 0.5e-12))
                else:                       #know c, guess l
                    self.guess_lclist.append(max(1/((2*pi*(1e6 + (40e6 - 1e6) * self.true_pkstarget[i]/(self.ptsnum - 1)))**2 * self.known_lclist[i]) - self.lxfmr2, 0.5e-6))

    def generate_pks_list(self, mxfmr, rlist, llist, clist):
        test = self.Z22(mxfmr, rlist, llist, clist)
        pkstest = find_peaks(np.real(test))[0]
        pkstestfit = []
        slicefit = []
        for pk in pkstest:
            x = [pk - 1, pk, pk + 1]
            y = np.real([test[x[0]], test[x[1]], test[x[2]]])
            mix, miy = threeptsfit(y, x)[1], threeptsfit(y, x)[2]
            pkstestfit.append(mix)
            slicefit.append(miy)
        return np.array(pkstestfit), np.array(slicefit)

    def prepare(self):
        def rmloss(p):
            guess_rlist = self.guess_rlist
            guess_rlist[-1] = p[0]
            if (np.min(guess_rlist) < 0):
                return np.array([1e6]*(self.ptsnum - self.thresh))
            return np.absolute(self.true_ttest - self.Z22(self.guess_mxfmr, guess_rlist, *self.translate_mask()))[self.thresh:]**2
        self.guess_rmloss = rmloss

        def rmloss2(p):
            guess_mxfmr = p[0]
            guess_rlist = p[1:]
            guess_rlist = np.append(guess_rlist, self.guess_rlist[-1])
            if (np.min(guess_rlist) < 0):
                return np.array([1e6]*self.thresh)
            return np.absolute(self.true_ttest - self.Z22(guess_mxfmr, guess_rlist, *self.translate_mask()))[:self.thresh]**2
        self.guess_rmloss2 = rmloss2

        def rmloss3(p):
            guess_rlist = self.guess_rlist
            guess_rlist[0] = p
            if (np.min(guess_rlist) < 0):
                return np.array([1e6]*self.thresh)
            return np.absolute(self.true_ttest - self.Z22(self.guess_mxfmr, guess_rlist, *self.translate_mask()))[:self.thresh]**2
        self.guess_rmloss3 = rmloss3

        def Clineloss(p):
            guess_Cline = p[0]
            guess_rlist = self.guess_rlist
            guess_llist, guess_clist = self.translate_mask()
            guess_clist = [*guess_clist[:self.loadnum], guess_Cline]
            if (p[0] < 0):
                return np.array([1e6]*(self.ptsnum - self.thresh))
            return np.absolute(self.true_ttest[self.thresh:] - self.Z22(self.guess_mxfmr, guess_rlist, guess_llist, guess_clist)[self.thresh:])**2
        self.guess_Clineloss = Clineloss
    
    def translate_mask(self): #translate from guesslc and known lc
        guess_llist = []
        guess_clist = []
        for i in range(self.loadnum + 1):
            if self.guess_mask[i] == 1:
                guess_llist.append(self.known_lclist[i])
                guess_clist.append(self.guess_lclist[i])
            else:
                guess_clist.append(self.known_lclist[i])
                guess_llist.append(self.guess_lclist[i])
        return guess_llist, guess_clist

    def translate_mask_p(self, guesslc, knownlc): #translate from customized guesslc and known lc
        guess_llist = []
        guess_clist = []
        for i in range(self.loadnum + 1):
            if self.guess_mask[i] == 1:
                guess_llist.append(knownlc[i])
                guess_clist.append(guesslc[i])
            else:
                guess_clist.append(knownlc[i])
                guess_llist.append(guesslc[i])
        return guess_llist, guess_clist

    def fitlc(self, smallC, jump0 = False):
        pkstot = self.true_pkstarget
        pksnow = self.generate_pks_list(self.guess_mxfmr, self.guess_rlist, *self.translate_mask())[0]
        currentguess = self.guess_lclist
        if jump0 == False:
            for i in range(5):
                currentguess[0] *= pksnow[0]/pkstot[0]
                pksnow = self.generate_pks_list(self.guess_mxfmr, self.guess_rlist, *self.translate_mask_p(currentguess, self.known_lclist))[0]
        for i in range(5):
            currentguess[1] *= pksnow[1]/pkstot[1]
            pksnow = self.generate_pks_list(self.guess_mxfmr, self.guess_rlist, *self.translate_mask_p(currentguess, self.known_lclist))[0]
            if jump0 == False:
                for j in range(5):
                    currentguess[0] *= pksnow[0]/pkstot[0]
                    pksnow = self.generate_pks_list(self.guess_mxfmr, self.guess_rlist, *self.translate_mask_p(currentguess, self.known_lclist))[0]
        for i in range(5):
            currentguess[2] *= pksnow[2]/pkstot[2]
            pksnow = self.generate_pks_list(self.guess_mxfmr, self.guess_rlist, *self.translate_mask_p(currentguess, self.known_lclist))[0]
            for j in range(5):
                if jump0 == False:
                    currentguess[0] *= pksnow[0]/pkstot[0]
                currentguess[1] *= pksnow[1]/pkstot[1]
                pksnow = self.generate_pks_list(self.guess_mxfmr, self.guess_rlist, *self.translate_mask_p(currentguess, self.known_lclist))[0]
        self.guess_lclist = currentguess

    def fitrm(self, loss):
        if loss == 3: # not consider rC and M
            self.guess_rlist[0] = 1/self.true_ttest[round(self.true_pkstarget[0])].real

    def fitC(self):
        Ccurrent = self.guess_lclist[-1]
        fittedC = leastsq(self.guess_Clineloss, [Ccurrent], maxfev=100)[0][0]
        self.guess_lclist[-1] = fittedC

    def solve(self, to_fitrc=False, jump0=False):
        if to_fitrc:
            self.fitlc(False, jump0)
            self.fitrm(3)
            return self.guess_lclist, self.guess_rlist
        else:
            self.fitlc(False, jump0)
            return self.guess_lclist

    def setparams(self, paramdic):
        self.lxfmr1 = paramdic["lxfmr1"]
        self.lxfmr2 = paramdic["lxfmr2"]
        self.guess_mask, self.known_lclist = paramdic["mask"]

class Csolver():
    def __init__(self):
        self.Rl = 0
        self.Csma = 0
        self.Lsma = 0
        self.Lxfmr1 = 0

    def calibrate(self, data):
        def calcS(p, x):
            Rl, Csma, Lsma, Lxfmr1 = p
            Z = 1/(1/(Rl + 2j * pi * x * Lxfmr1) + 1/(2j * pi * x * Lsma + 1/(2j * pi * x * Csma)))
            S = (Z - 50)/(Z + 50)
            return S

        def caliloss(p, y, x):
            S = calcS(p, x)
            return np.absolute(y - S)**2

        arx = np.linspace(1e6, 40e6, 101)
        arcomplex = data
        aromega = 2 * pi * arx
        ptsnum = 101
        arZ = Zfunc(arcomplex)
        arY = 1/arZ
        Yreal = np.real(arY)
        Zreal = 1/Yreal
        ix = (0, int(0.05*ptsnum), int(0.1*ptsnum))
        ixx = (arx[ix[0]], arx[ix[1]], arx[ix[2]])
        ixy = (Zreal[ix[0]], Zreal[ix[1]], Zreal[ix[2]])
        k = diff_2(ixy, ixx)
        guess_R_l = Zreal[ix[0]] - k * arx[0]**2
        guess_L_xfmr1 = sqrt(k * guess_R_l/(4 * pi**2))
        Yimag = np.imag(arY - 1/(guess_R_l + 1j * guess_L_xfmr1 * aromega))
        idx = int(0.5*ptsnum)
        guess_C_SMA = np.mean(Yimag[:idx]/aromega[:idx])
        guess_L_SMA = (1/(aromega * guess_C_SMA) - 1/Yimag)[-1]/aromega[-1]
        guessparam = (guess_R_l, guess_C_SMA, guess_L_SMA, guess_L_xfmr1)
        fittedparam = leastsq(caliloss, guessparam, args=(arcomplex, arx), maxfev=100)
        self.Rl, self.Csma, self.Lsma, self.Lxfmr1 = fittedparam[0]

    def unpack(self, data):
        '''
        Z2 = Z12^2 / Z22
        即副电路的阻抗
        VNA上的值和副电路是一一对应的
        和calibrate状态的差距可以看作副电路的值
        '''
        arx = np.linspace(1e6, 40e6, 101)
        arcomplex = data
        Z = Zfunc(arcomplex)
        aromega = 2 * pi * arx
        Zcma = 1j * aromega * self.Lsma + 1/(1j * aromega * self.Csma)
        Z1 = self.Rl + 1j * aromega * self.Lxfmr1
        Z2 = 1/(1/Z - 1/Zcma) - Z1
        return Z2