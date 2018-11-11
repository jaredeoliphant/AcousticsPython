import numpy as np
from acousticsFunctions import binfileload
from spectra import autospec, fractionalOctave
import sys
import matplotlib.pyplot as plt

path = sys.path[0]+'/ReverbFiles'
print (path)

N = 6195200
x1 = np.zeros((N,6))
x2 = x1
for i in range(6):
    x1[:,i] = binfileload(path,'ID',1,i,N)
    x2[:,i] = binfileload(path,'ID',2,i,N)

fs = 102.4e3
dt = 1/fs
T = 60.5
t = np.arange(0,T,dt)
pref = 2e-5
ns = 2**15
unitflag = 0
prop_distance = ns/fs*343
print("In 1 block we travel %.2f meters" %prop_distance)
print("Frequency resolution is %.0f Hz" %(fs/ns))


Gxx1 = np.zeros((int(ns/2),6))
Gxx2 = Gxx1
spec1 = np.zeros((31,6))
spec2 = spec1
fig1, ax1 = plt.subplots()
for i in range(6):
    Gxx1[:,i], f, OASPL = autospec(x1[:,i], fs, ns, N, unitflag)
    Gxx2[:,i], f, OASPL = autospec(x2[:,i], fs, ns, N, unitflag)
    spec1[:,i], fc = fractionalOctave(f,Gxx1[:,i])
    spec2[:,i], fc = fractionalOctave(f,Gxx2[:,i])
    ax1.semilogx(fc,10*np.log10(spec1[:,i]/pref**2))


ax1.set_xlabel("Frequecy 1/3 octave bands (Hz)")
ax1.set_ylabel("SPL (dB re 20$\mu$Pa)")
ax1.set_xlim((100,10e3))

fc = fc[7:28]
spec1 = spec1[7:28,:]
spec2 = spec2[7:28,:]

Lp_bar1 = np.zeros((len(spec1),1))
Lp_bar2 = Lp_bar1
for i in range(len(spec1)):
    Lp_bar1[i] = 10*np.log10(np.sum(spec1[i,:]/pref**2)/6)
    Lp_bar2[i] = 10*np.log10(np.sum(spec2[i,:]/pref**2)/6)

fig2, ax2 = plt.subplots()
ax2.semilogx(fc,Lp_bar1)
ax2.semilogx(fc,Lp_bar2)




# T60 numbers from Travis for Large chamber
f = np.array([100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000],dtype=float)
T60 = np.array([8.9975,8.663333333,7.614166667,7.469166667,7.964166667,8.323333333,8.4825,8.195,7.944166667, \
7.798333333,7.245,6.283333333,5.4575,4.64,3.7775,3.179166667,2.564166667,1.899166667,1.375833333,1.079166667,0.6941666667],dtype=float)

fig3, ax3 = plt.subplots()
ax3.semilogx(f,T60)


# Metetoriological
T = 30.0   # Temperature in Celsius
B = 101340.0  # Barometric Pressure in Pa
B0 = 1.013e5  # reference Pressure Pa

# speed of sound at temperature T
c = 20.05*np.sqrt(273+T)  # m/s

# Room properties
V = 4.96 * 5.89 * 6.98 #  volume of the room (m^3)
A = 55.26*V/c/T60 # equivalent absorbption area of the room as function of freq (m^2)
A0 = 1.0  # m^2
S = 2*(4.96*5.89) + 2*(5.89*6.98) + 2*(6.98*4.96) # total surface area of the room (m^2)

## check absorption requirements (5.3)
# fprintf('The number of frequecies that meet the absorption requirements\nis %d/%d. Trev > V/S = %.2f\n',nnz(T60 > V/S), length(T60),V/S)

"""
%% dmin
C1 = 0.08;
C1 = 0.16; % to minimize near-field bias error
dmin = min(C1*sqrt(V./T60)) % minimum distance b/t source and microphone (m)
% the microphones shall be more that 1.0 m from a wall
% the min distance between mics is half the wavelength of the lowest freq.
% or interest (100 Hz = 3.4/2 = 1.7 meters)
% with 1.1 m spacing I can go down to 156 Hz
"""
########Final Equation!
# sound power level of the source as a function of frequency
Lw1 = Lp_bar1.transpose()[0] + 10*np.log(A/A0) + 4.34*A/S + 10*np.log10(1+S*c/8/V/f) - 25*np.log10(427*np.sqrt(273.0/(273+T))*B/B0/400.0) - 6
Lw2 = Lp_bar2.transpose()[0] + 10*np.log(A/A0) + 4.34*A/S + 10*np.log10(1+S*c/8/V/f) - 25*np.log10(427*np.sqrt(273.0/(273+T))*B/B0/400.0) - 6

fig4, ax4 = plt.subplots()
ax4.plot(Lw1)  #,'color',[.75 .6 0],'linewidth',8,'marker','none')
ax4.plot(Lw2)  #,'color',[0 0 .75],'linewidth',8,'marker','none')
# ax4.set_xscale('log')
# xlim([90 11e3])
# % ylim([0 10])
# set(gca,'xminortick','off')
# set(gca,'tickdir','out')
# set(gca,'xtick',f)
# set(gca,'xticklabels',{'','125','','','250','','','500',...
#     '','','1000','','','2000','','','4000','','','8000',''})
# ylabel('L_w (dB re 1pW)')
# xlabel('1/3rd-octave freq. band (Hz)')
# grid off
# box off
# legend('1st mic positions','2nd mic positions','location','northwest')
# print('sound_power_2times.tif','-dtiff','-r300')
# %% Standard deviation (dB sense) (calculate for each frequency band)
# NM = 6; % number of microphones

# for i = 1:length(Lp_bar1)  % loop through each freq.
#     Lp_mics = 10*log10(spec1(i,:)/pref^2);  % SPL of all the mics at a given freq
#     Lpm = Lp_bar1(i); % the arithmetic mean of the SPL for 6 microphones
#     sM(i) = sqrt(sum((Lp_mics-Lpm).^2/(NM - 1)))
#     % if sM < 1.5 for all freq. bands you are good!
#     sM < 1.5
# end



plt.show()
