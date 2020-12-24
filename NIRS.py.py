import re
import math
import sys

import numpy as np
import scipy
from scipy.signal import butter, lfilter, filtfilt, freqz, firwin #Band-Pass

import matplotlib.pyplot as plt

import tkinter as tk # GUI
from tkinter import ttk
from tkinter import messagebox as mb
from tkinter import simpledialog as sd # Input string and value
from tkinter import filedialog as fd # Read string and value, and open file path

import PySide2 as ps2  # <<< LGPL
from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *

class PreProcessing():
    def __init__(self):
        pass

    
    #---------
    # Parameter extraction
    #---------
    def references(self):
        '''Select the file, and convert the data to a list, and finally draw the picture'''
        # Open window
        loadFile = tk.Tk()
        loadFile.withdraw() # Hide tk window
        self.fileman = fd.askopenfilename() # Unexplained delay
        self.fileman = self.fileman[self.fileman.find(":")-1:self.fileman.find(".")] # File path, go to extension
        intro = open(self.fileman +'.hdr', 'r')
        info = intro.read()
        intro.close
        # Read in all parameter settings
        GeneralInfo = info[info.find("GeneralInfo") +12:info.find("ImagingParameters")-1]
        ImagingParameters = info[info.find("ImagingParameters")+18:info.find("Paradigm")-1]
        length1 = int(ImagingParameters[ImagingParameters.find("Wavelengths")+13:ImagingParameters.find("Wavelengths")+16])
        length2 = int(ImagingParameters[ImagingParameters.find("Wavelengths")+17:ImagingParameters.find("TrigIns")-2])
        self.Length = [length1, length2]
        self.Fs = float(ImagingParameters[ImagingParameters.find("SamplingRate")+13:ImagingParameters.find("Mod Amp")])
        Paradigm = info[info.find("Paradigm")+9:info.find("ExperimentNotes")-1]
        ExperimentNotes = info[info.find("ExperimentNotes")+16:info.find("GainSettings")-1]
        # Gains matrix
        GainSettings = info[info.find("GainSettings")+22:info.find("Markers")-5]
        self.GainSettings = list(map(int,re.findall(r"\d|-\d", GainSettings)))
        # triggers, Conditions (in frames)
        Markers = info[info.find("Markers")+19:info.find("DataStructure")-5]
        # If you want to generate str matrix, Markers = Markers.strip("\n").strip("\t")
        Markers = Markers.split()
        # Split out frames、triggers、seconds
        self.seconds = []
        self.triggers = []
        self.frames = []
        for i in range(len(Markers)):
            if (i+1) % 3 == 0:
                self.frames.append(Markers[i])
            elif (i+2) % 3 == 0:
                self.triggers.append(Markers[i])
            else:
                self.seconds.append(Markers[i])
        # Data structure
        DataStructure = info[info.find("DataStructure")+14:info.find("DarkNoise")-5]
        S_D_Key = re.split(r":|,", DataStructure[DataStructure.find("S-D-Key")+9:DataStructure.find("S-D-Mask")-3])  # 多重切割
        self.S_D_Keys = []
        S_D_Values = []
        for i in range(len(S_D_Key)):
            if i % 2 == 0:
                self.S_D_Keys.append('S'+S_D_Key[i].replace('-', "_"))
            else:
                S_D_Values.append('S'+S_D_Key[i].replace('-', "_"))  # Is it fixed as a serial number?
        self.S_D_Mask = DataStructure[DataStructure.find("S-D-Mask")+11:len(DataStructure)].split()
        # Dark Noise
        DarkNoise = info[info.find("DarkNoise")+10:info.find("ChannelsDistance")-1]
        # Channels Distance
        ChannelsDistance = info[info.find("ChannelsDistance")+13:len(info)].split() # Is the length of the available Channel (d)
        # Read data
        wl1 = open(self.fileman+'.wl1', 'r')
        wl2 = open(self.fileman+'.wl2', 'r')
        data1 = wl1.read().split()  # Total length = Channels(1.) x frames()
        data2 = wl2.read().split()
        self.max_frames = int(len(data1) / len(self.S_D_Mask))  # Is not an integer foolproof?
        self.total_channels = int(len(data1) / self.max_frames)
        wl1.close()
        wl2.close()
        # Dimension conversion
        self.xLabel = []
        self.raw_wl1 = []
        self.raw_wl2 = []
        collect_1_DataPoint = []
        collect_2_DataPoint = []
        for i in range(self.total_channels):  # 512
            for j in range(self.max_frames): # 868
                if i == 0:
                    self.xLabel.append(j+1) # Generate x axis
                collect_1_DataPoint.append(data1[i + self.total_channels*j])
                collect_2_DataPoint.append(data2[i + self.total_channels*j])
                if j == self.max_frames - 1:
                    self.raw_wl1.append(list(map(float,collect_1_DataPoint))) # As a result, str is converted to float
                    self.raw_wl2.append(list(map(float,collect_2_DataPoint)))
                    collect_1_DataPoint = []
                    collect_2_DataPoint = []
    #-------------
    # Global line graph
    #-------------
    def raw_lineplot(self):
        self.raw = {
            "self.raw_wl1": ["self.Length[0]",1,"self.filtered_wl1"], 
            "self.raw_wl2": ["self.Length[1]",2,"self.filtered_wl2"]
            } # {Original line graph: [wavelength, location of image sub-version, filtered line graph]}
        for key,value in self.raw.items():
            for i in range(len(self.S_D_Mask)):  # 0 1
                if self.S_D_Mask[i] == "1":
                    fig = plt.subplot(2,1,value[1])
                    plt.plot(self.xLabel, eval(key)[i])#, label=self.S_D_Keys[i]) 
                    # raw data
            fig.title.set_text("Available Channels ("+ str(eval(value[0])) +"nm)")
            plt.grid(True)
            plt.ylabel('Amplitude(V)')
        plt.xlabel('frames('+str(round(float(self.Fs),2))+'Hz)')
        plt.show()
    #----------------------------------
    # Discontinuities, Spike Artifacts
    #----------------------------------
    def head_displace(self):
        '''Manually mark, clear all channels, start processing from the largest spike, both sides to the middle, *triggers in the range are not counted*'''
        reject_band = 0
        reject_area = []
        while reject_band != "z" or reject_band == None:
            temp0 = 0.00
            temp1 = 0.00
            temp2 = 0.00
            f0 = 0
            f1 = 0
            f2 = 0
            try:
                reject_band = input("Enter the start and end points (frames, as integers) of the recognized head movement, separated by a slash, such as "28/95", and end input ：")
                reject1 = int(reject_band[0:reject_band.find("/")])
                reject2 = int(reject_band[reject_band.find("/")+1:len(reject_band)])
                reject_area.append(reject1) # <<<< Need sorting algorithm
                reject_area.append(reject2)
                for raw in [self.raw_wl1,self.raw_wl2]:
                    for i in range(len(self.S_D_Mask)):
                        if reject1 < 2: # Foolproof
                            reject1 = 2
                        if reject2 > self.max_frames:
                            reject2 = self.max_frames
                        if reject1 > reject2:
                            reject1, reject2 = reject2, reject1 # Value swap
                        if reject2 > len(raw[i]):
                            reject1 = len(raw[i])
                        if self.S_D_Mask[i] == "1":
                            for j in range(0,reject1-1): # before
                                temp1 += raw[i][j]
                                f1 += 1
                            for j in range(reject1-1,reject2-1): # in
                                temp0 += raw[i][j]
                                f0 += 1
                            for j in range(reject2-1,len(raw[i])): # rear
                                temp2 += raw[i][j]
                                f2 += 1
                            rangeX = reject2 - reject1
                            rangeY = temp2/f2 - temp1/f1
                            slope = rangeY/rangeX # Consequential slope
                            if f1 >= f2 : # The rear section is adapted to the front section
                                for j in range(reject2-1, len(raw[i])): 
                                    raw[i][j] = raw[i][j] * (temp1/f1)/(temp2/f2) # <<<< Message compression
                                for j in range(reject1-1, reject2-1): # Reject the correction in the domain <haven't found a good way
                                    raw[i][j] = temp1/f1
                            elif f2 > f1: # Front section and back section
                                for j in range(0, reject2-1): 
                                    raw[i][j] = raw[i][j] * (temp2/f2)/(temp1/f1)
                                for j in range(reject1-1, reject2-1): # Reject corrections in the domain
                                    raw[i][j] =temp2/f2
                            temp0 = 0.00
                            temp1 = 0.00
                            temp2 = 0.00
                            f0 = 0
                            f1 = 0
                            f2 = 0
                PreProcessing.raw_lineplot(self)
            except:
                print("try again, or ?")

    #----------------------------------
    # Coeffiecient of Variation(CV)
    # According to the standard Good or Bad Channels
    #----------------------------------
    
    # CV(Coefficient of Variation) = S / Xbar
   
    def cov(self):
        temp1 = 0.00 # Channel value summation
        temp2 = 0.00
        temp1sq = 0.00 # Channel value sum of squares
        temp2sq = 0.00
        wl1_avg = [] # Generate the average and standard deviation of each Channel
        wl1_delta = []
        wl2_avg = []
        wl2_delta = []
        self.cv1 = [] # Coefficient of Variation
        self.cv2 = []

       # First calculate the overall mean standard deviation (parent group)
        for i in range(len(self.S_D_Mask)):  # 0 1
            if self.S_D_Mask[i] == "1":
                for j in range(len(self.raw_wl1[i])):
                    w1 = float(self.raw_wl1[i][j])
                    w2 = float(self.raw_wl2[i][j])
                    temp1 += w1
                    temp2 += w2
                    temp1sq += w1**2
                    temp2sq += w2**2
                m1 = temp1/self.max_frames
                s1 = pow((temp1sq-temp1**2/self.max_frames)/self.max_frames, .5)
                m2 = temp2/self.max_frames
                s2 = pow((temp2sq-temp2**2/self.max_frames)/self.max_frames, .5)
                wl1_avg.append(m1)
                wl2_avg.append(m2)
                wl1_delta.append(s1)
                wl2_delta.append(s2)
                self.cv1.append( 100*s1/m1) # %
                self.cv2.append( 100*s2/m2)
                temp1 = 0.00
                temp2 = 0.00
                temp1sq = 0.00
                temp2sq = 0.00
        #return cv1, cv2
    # gain_Setting = int(input("Gain standard："))
    # channel_CV = float(input("Channel's C.V. (%)："))

    #-----------------------
    # Bandpass Filter
    #-----------------------
    def butter_bandpass(self, data, lowcut, highcut, fs, order=5):
        '''
        
        Butterworth filter
        . fs Sampling rate
        The order of the filter, the larger the order, the better the effect, but the amount of calculation also increases
        The cutoff frequency of Wn normalization is between 0-1, and the highest frequency that can be processed is fs/2
        Nirslab > LPF = .01, HPF = .2
        Lu, Chia-Feng > LPF = .01, HPF = .1
        '''
        # The Nyquist rate of the signal.
        #nyq = .5 * fs 
        #low = lowcut / nyq # .01
        #high = highcut / nyq # .2

        b, a = butter(order, [lowcut,highcut], btype='band')
        #b, a = butter(order, lowcut, btype='low')
        #b, a = butter(order, highcut, btype='high')

        #w, h = freqz(b, a, worN=2000) # Add
        #y = lfilter(b, a, data) 
        y = filtfilt(b, a, data) # < Forward and backward twice
        return y

    def filterLine(self):
        '''
        From the original line graph> filtered line graph (presented)> blood oxygen concentration line graph (A)
            raw_        filtered_    deltaA_
        
        lf = .01, hf = .2 in nirsLAB
        But the actual test results of lf = .005 and hf = .12 are relatively close 
        # not divided by nyq?
        '''
        temp1 = 0.00
        temp2 = 0.00
        self.A1 = []
        self.A2 = []
        self.filtered_wl1 = self.raw_wl1
        self.filtered_wl2 = self.raw_wl2
        for key,value in self.raw.items():
            fig = plt.figure() # Empty the canvas
            fig = plt.subplot(111)
            for i in range(len(self.S_D_Mask)):  # 0 1
                if self.S_D_Mask[i] == "1": 
                    fig = plt.subplot(111)
                    if value[1] == 1:
                        for j in range(len(self.raw_wl1[i])):
                            temp1 += self.raw_wl1[i][j]
                        temp1 = temp1/len(eval(key)[i]) # Average
                        self.filtered_wl1[i] = PreProcessing.butter_bandpass(self,self.filtered_wl1[i], .005, .12, fs=self.Fs) #Band-Pass
                        self.filtered_wl1[i] = self.filtered_wl1[i] + temp1 # For presentation
                        plt.plot(self.xLabel, self.filtered_wl1[i], label=self.S_D_Keys[i]) # filtered data
                        temp1 = 0.00
                    elif value[1] == 2:
                        for j in range(len(self.raw_wl2[i])):
                            temp2 += self.raw_wl2[i][j]
                        temp2 = temp2/len(self.raw_wl2[i]) # Average
                        self.filtered_wl2[i] = PreProcessing.butter_bandpass(self,self.filtered_wl2[i], .005, .12, fs=self.Fs) #Band-Pass
                        self.filtered_wl2[i] = self.filtered_wl2[i] + temp2 # For presentation
                        plt.plot(self.xLabel, self.filtered_wl2[i], label=self.S_D_Keys[i]) # filtered data
                        temp2 = 0.00
            #plt.ylim(-0.003, 0.003)
            fig.title.set_text("Filtered Data("+str(eval(value[0]))+"nm)")
            plt.grid(True)
            plt.xlabel('frames')
            plt.ylabel('Amplitude(V)')
            plt.legend(loc='best') # show label
            plt.pause(.1)
            plt.show()
    #-------------------------
    # Beer-Lambert Law (mBLL)
    #-------------------------    
    
    # Molar Extinction Coefficients [wavelength, for oxyH, for deoxyH]
    
    def concentration(self): # input filtered data
       
        wave_Lengths = 2 # 760 850
        baseline_inframes = [1, self.max_frames]
        totHb = 75 # uM
        mov2Sat = 70 # %

        pathLength = 3 # Measuring distance, optical path(cm), the depth is equivalent to 1/4 of the distance
        wl1_DPF = 6.4 # Differential Pathlength Factor (DPF)
        wl2_DPF = 5.75

        dp1 = wl1_DPF * pathLength # differential pathlength(DP)
        dp2 = wl2_DPF * pathLength

        sq = 10**(-3)
        εHbO_wl1 = 1486.5865*sq # unit
        εHbO_wl2 = 2526.3910*sq
        εHbR_wl1 = 3843.707*sq
        εHbR_wl2 = 1798.643*sq
        
        #----------------------------------------------
        # Light absorbance (A) = absorption coefficient (α) * optical path length (l) * concentration (c)
        #----------------------------------------------
   
        #--------
        # I to A 
        #--------

        temp_A1 = []
        temp_A2 = []
        self.deltaA_wl1 = []
        self.deltaA_wl2 = []
        for i in range(len(self.S_D_Mask)):  # 0 1
            if self.S_D_Mask[i] == "1":
                for j in range(self.max_frames):
                    temp_A1.append(math.log(np.mean(self.filtered_wl1[i]) / self.filtered_wl1[i][j],math.e))
                    temp_A2.append(math.log(np.mean(self.filtered_wl2[i]) / self.filtered_wl2[i][j],math.e))
                self.deltaA_wl1.append(temp_A1)
                self.deltaA_wl2.append(temp_A2)
                temp_A1 = []
                temp_A2 = []
            else:
                self.deltaA_wl1.append([])
                self.deltaA_wl2.append([])
        #---------------
        # A to HbO & Hb
        #---------------
        b = float(pathLength*(εHbR_wl1*εHbO_wl2 - εHbO_wl1*εHbR_wl2))
        temp_HbO = []
        temp_HbR = []
        self.HbO = []
        self.HbR = []
        for i in range(len(self.S_D_Mask)):  # 0 1
            if self.S_D_Mask[i] == "1":
                for j in range(self.max_frames):
                    temp_HbO.append((εHbR_wl1*self.deltaA_wl2[i][j]/wl2_DPF - εHbR_wl2*self.deltaA_wl1[i][j]/wl1_DPF)/b)
                    temp_HbR.append((εHbO_wl2*self.deltaA_wl1[i][j]/wl1_DPF - εHbO_wl1*self.deltaA_wl2[i][j]/wl2_DPF)/b)
                self.HbO.append(temp_HbO)
                self.HbR.append(temp_HbR)
                temp_HbO = []
                temp_HbR = []
            else:
                self.HbO.append([])
                self.HbR.append([])

        return self.HbO, self.HbR

    def hemoglobin(self):
        PreProcessing.concentration(self)
        fig = plt.figure() # Empty canvas
        fig = plt.subplot(111)
        for i in range(0,5):  #len(self.S_D_Mask)):  # 0 1
            if self.S_D_Mask[i] == "1": 
                plt.plot(self.xLabel, self.HbO[i], label=self.S_D_Keys[i]) # filtered data
                plt.plot(self.xLabel, self.HbR[i], label=self.S_D_Keys[i]) # filtered data
        #plt.ylim(-0.003, 0.003)
        plt.grid(True)
        #plt.xlabel('frames')
        #plt.ylabel('Amplitude(V)')
        plt.legend(loc='best') # show label
        plt.pause(.1)

#---------------------
# GUI >>> Initial input interface
#---------------------
pp = PreProcessing()
workshop = tk.Tk()

workshop.title("NiRx data analysis")
workshop.resizable(0,0) # Lock window size

def clickOK():
    mb.showinfo("Analysis begins", "Ready to start! \n Let's get out of hand!")

def mb_Option():
    mb.showinfo()
    mb.showerror()
    mb.askyesno()
    mb.askokcancel()

def cutpoint0():
    sd.askinteger("Cut Point", "Please enter the starting cut point(frames, such as-10)：")

def cutpoint1():
    sd.askinteger("Cut Point", "Please enter the end cut point(frames, Such as 80)：")

# Introduction
preProcessing = tk.Label(workshop,bg="red",text="PreProcessing",fg="white", width=10,height=2,borderwidth=10).grid(column=0,row=0)
level1 = tk.Label(workshop,bg="green",text="LEVEL 1.",fg="white",width=10,height=2,borderwidth=10).grid(column=1,row=0)
level2 = tk.Label(workshop,bg="blue",text="LEVEL 2.",fg="white",width=10,height=2,borderwidth=10).grid(column=2,row=0)

buttons_func = {
    ## PreProcessing
    # Data & Conditions
    "load_Data":     ["Get Parameters", 0, 1, pp.references],
    "set_Markers":   ["Image-Signal Graphic", 0, 2, pp.raw_lineplot],
    # Data Preprocessing 
    "truncate":      ["Move", 0, 3, pp.head_displace],
    "check_Quality": ["Covariation analysis", 0, 4, clickOK],
    "apply_Filter":  ["Draw a filter map", 0, 5, pp.filterLine],
    # Hemodynamic Staties
    "set_Parameters":["Try result", 0, 6, pp.hemoglobin],
    "compute":       ["Context archive", 0, 7, clickOK],
    ## First level
    "first_lv":      ["Individual comparison", 1, 1, clickOK],
    ## Second level
    "second_lv":     ["Experiment comparison", 2, 1, clickOK],
    }

for key in buttons_func:
    key = ttk.Button(workshop, 
        text=buttons_func[key][0], 
        command=buttons_func[key][3]).grid(column=buttons_func[key][1],row=buttons_func[key][2],sticky="ew") 

workshop.mainloop() # Call and maintain the window