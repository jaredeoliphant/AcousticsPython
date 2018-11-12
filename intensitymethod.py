import numpy as np
import sys
from acousticsFunctions import binfileload, weighting
from spectra import autospec,crossspec, fractionalOctave
import matplotlib.pyplot as plt

# path to the files of interest
side = 2
path = sys.path[0]+"/IntensityFiles/Side"+str(side)
print ("path to files: ",path)

# recording information from log file
fs = 50000.0
dt = 1/fs
T = 10.5
N = int(fs*T)
t = np.arange(0,T,dt)


# one of the sides only has 79 recordings (id's)
idnums = 81
if side == 2:
    print("Only 79 files for this side")
    idnums = 79

# initialize the 2d arrays for the two microphones
x = np.zeros((N,idnums))
y = np.zeros((N,idnums))

print("loading in the data...")
for i in range(idnums):   # looping through the different ID numbers (81 per side)
    # each "column" of y and x is a different ID
    y[:,i] = binfileload(path,'ID',i+1,0,N)  # farther mic to the source 
    x[:,i] = binfileload(path,'ID',i+1,1,N)  # closer mic to the source


print("2 arrays built with shape: ", np.shape(x))


# pressure and power/intensity references
pref = 2e-5
iref = 1e-12

# intensity calculation parameters
ns = 2**15   # samples per block
rho = 1.2    # denisty of the air (find based on temp, pressure, RH)
deltax = .0254  # spacing between microphones
Area1 = 0.15*0.15   # area of one measurement (15 cm distance traveled)
Area2 = (1.2+0.15)**2  # total area of one side

print("Areas match? ",Area2 - 81*Area1 < .0001) # Area2 should be 81*Area1


fig1,ax1 = plt.subplots(figsize=(10,10)) 
Isum = 0
Iavg_over_freq = np.zeros((idnums,1))
# loop through all the IDs and sum up intensity as we go
print("Calculating Intensity...")
for i in range(idnums):
    Gxy,f = crossspec(x[:,i],y[:,i],fs,ns,N) # crossspec for each ID (column)
    f = f[1:]     # cut out zero Hz
    Gxy = Gxy[1:]
    Intensity = np.imag(Gxy) / 2.0 / np.pi / f / rho / deltax    # equation for intensity
    Iavg_over_freq[i] = np.mean(np.log10(np.abs(Intensity)/iref))
    Isum += Intensity * Area1     # sum up over total area
    ax1.semilogx(f,10*np.log10(np.abs(Intensity)/iref))   # plot each id seperate
    
ax1.set_xlabel("Frequency (Hz)")
ax1.set_ylabel("Intensity (dB re 1pW/m$^2$)")
ax1.set_title("Intensity from each recording")

Iavg = Isum / Area2 #surface sound intensity for one side (The areas should do nothing in this case because each measurement square is identical)
print("Iavg over the surface is ",Iavg)
fig2, ax2 = plt.subplots(figsize=(10,10))
ax2.semilogx(f,10*np.log10(np.abs(Iavg)/iref))
ax2.set_xlabel("Frequency (Hz)")
ax2.set_ylabel("Intensity (dB re 1pW/m$^2$)")
ax2.set_title("Average Intensity for side "+str(side))



## make a plot showing OASPL vs IDnum
fig3, ax3 = plt.subplots(figsize=(10,10))
ax3.plot(Iavg_over_freq)
ax3.set_title("average intensity versus id number")






## convert to a single value to be reported as the A-weighted sound power level
ff = fractionalOctave(np.arange(2,4),np.arange(3,5),flims=[100,20e3],width=3)[1]  # get the 1/3 octave band freq. out
Gain = weighting(ff,type='A')[1]  # only save the second output in this case
#Overall Sound power level
#Lw_overall = 10*np.log10(np.sum(10**(.1*(Lw+Gain))))   # where C is the A-weighting constant  
# print()
# print("The A-weighted overall sound power level is: ",Lw_overall)
print()
plt.show()
