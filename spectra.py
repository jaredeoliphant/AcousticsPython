"""
Module 'spectra' contains autospec, crossspec, and fractionalOctave



This program calulates the autospectral density or autospectrum and the OASPL of a signal.
Hanning windowing is used, with 50% overlap. Per Bendat and Piersol, Gxx 
is scaled by the mean-square value of the window to recover the correct OASPL.
call [Gxx,f,OASPL] = autospec(x,fs,ns,N);
Outputs: 
Gxx = Single-sided autospectrum or autospectral density, depending on unitflag
f = frequency array for plotting
OASPL = Overall sound pressure level
Inputs:
x = time series data.
fs = sampling frequency
ns = number of samples per block.  Default is 2^15 if not specified.
N = total number of samples.  If N is not an integer multiple of ns, 
the samples less than ns in the last block are discarded.  Default   
is nearest lower power of 2 if not specified.
unitflag = 1 for autospectrum, 0 for autospectral density.  Default is
autospectral density
Authors: Kent Gee, Alan Wall, and Brent Reichman
translated to python by Jared Oliphant
"""

import numpy as np 
import math

def autospec(x,fs,ns,N,unitflag=0):

    # coerce inputs
    ns = int(ns)

    # frequency array
    f = (fs/ns)*np.arange(0,ns/2.)
    df = f[1]

    # enforce zero mean
    x -= np.mean(x)

    ww = np.hanning(ns)
    W = float(np.mean(ww**2))

    numBlocks = int(math.floor(2*N/ns-1))

    blockMat = np.tile(np.matrix(np.arange(0,ns,dtype=int)),[numBlocks,1]) \
     + np.tile(np.matrix(np.arange(0,numBlocks,dtype=int)).T*ns/2,[1,ns])
    blockMat = np.matrix(blockMat, dtype=int)

    blocks = np.multiply(np.tile(ww,[numBlocks,1]),x[blockMat])
    X = np.fft.fft(blocks)
    
    Xss = X[:,0:int(ns/2)]

    Scale = 2/float(ns)/fs/W
    Gxx = Scale*np.mean(np.conjugate(Xss)*Xss,axis=0)

    Gxx = Gxx*df**unitflag
    Gxx = np.real(Gxx)

    if unitflag == 0:
        OASPL = 20*np.log10(np.sqrt(np.sum(Gxx*df))/2e-5)
    else:
        OASPL = 20*np.log10(np.sqrt(np.sum(Gxx))/2e-5)

    return [Gxx,f,OASPL]


















'''
This program calulates the crossspectral density or spectrum of signals x and y.
Hanning windowing is used, with 50% overlap. Per Bendat and Piersol, Section 11.6.3, Gxy 
is scaled by the mean-square value of the window for overall amplitude
scaling purposes.
call [Gxy,f] = crossspec(x,y,fs,ns,N,unitflag);
Outputs: 
Gxy = Single-sided cross spectrum or cross spectral density, depending on unitflag
f = frequency array for plotting
Inputs:
x,y = time series data
fs = sampling frequency
ns = number of samples per block.  Default is 2^15 if not specified.
N = total number of samples.  If N is not an integer multiple of ns, 
the samples less than ns in the last block are discarded.  Default   
is nearest lower power of 2 if not specified.
unitflag = 1 for autospectrum, 0 for autospectral density.  Default is
autospectral density
Authors: Kent Gee and Alan Wall; Translated to python by Jared Oliphant
'''
import numpy as np 
import math

def crossspec(x,y,fs,ns=2**15,N=-1,unitflag=0):

    if N == -1:
        N = 2**math.floor(np.log2(len(x)))

    
    # frequency array
    f = (fs/ns)*np.arange(0,ns/2.)
    df = f[1]

    # enforce zero mean
    x -= np.mean(x)
    y -= np.mean(y)

    ww = np.hanning(ns)
    W = float(np.mean(ww**2))

    numBlocks = math.floor(2*N/ns-1)

    blockMat = np.tile(np.matrix(np.arange(0,ns,dtype=int)),[numBlocks,1]) \
     + np.tile(np.matrix(np.arange(0,numBlocks,dtype=int)).T*ns/2,[1,ns])


    blocksx = np.multiply(np.tile(ww,[numBlocks,1]),x[blockMat])
    blocksy = np.multiply(np.tile(ww,[numBlocks,1]),y[blockMat])
    X = np.fft.fft(blocksx)
    Y = np.fft.fft(blocksy)

    Xss = X[:,0:ns/2]
    Yss = Y[:,0:ns/2]

    Scale = 2/float(ns)/fs/W
    Gxy = Scale*np.mean(np.multiply(np.conjugate(Xss),Yss),axis=0)

    Gxy = Gxy*df**unitflag
    
    return [Gxy,f]


















'''
[spec,fc]=FractionalOctave(f,Gxx,flims,width)
Performs a frequency-domain, fractional-octave analysis using filter
masks using the ANSI 2004 standard. Filter masks calculated using exact
center frequencies (referenced to 1 kHz), whereas preferred frequencies
are returned.
Inputs:   f - frequency array (Hz)
Gxx - autospectral density in Engineering Units^2/Hz
flims - [flow, fhigh], desired range of low and high frequency
fractional-octave bands between 1e-2 and 1e6 Hz.
Default is [20,20000];  User should ensure the lowest
frequency selected is a preferred center frequency for
the selected bandwidth.
width - fractional octave bandwidth, 1/width. Options are
1,3,6,12,and 24. Default is width=3;
Outputs:  fc, preferred band center frequencies
spec, octave band spectra (Eng Units^2)
Authors: Kent Gee
'''
import numpy as np

def fractionalOctave(f,Gxx,flims=[2e1,2e4],width=3):

    fcsub = np.array([1,1.03,1.06,1.09,1.12,1.15,1.18,1.22,1.25,1.28,\
    1.32,1.36,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2,2.06,\
    2.12,2.18,2.24,2.3,2.36,2.43,2.5,2.58,2.65,2.72,2.8,2.9,3,3.07,3.15,3.25,\
    3.35,3.45,3.55,3.65,3.75,3.87,4,4.12,4.25,4.37,4.5,4.62,4.75,4.87,5,5.15,5.3,\
    5.45,5.6,5.8,6,6.15,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.75,8,8.25,8.5,8.75,9,9.25,9.5,9.75])

    fc= np.concatenate((fcsub*1e-2,fcsub*1e-1,fcsub,fcsub*1e1,fcsub*1e2,fcsub*1e3,fcsub*1e4,\
    fcsub*1e5))
    fc = np.append(fc,1e6)


    allowwidths = [1,3,6,12,24]
    if width not in allowwidths:
        print('bad width')
    

    n = np.arange(0,len(fcsub))
    fcsubexact = 1000*2**(n/24.)
    fcexact= np.concatenate((fcsubexact*1e-5,fcsubexact*1e-4,fcsubexact*1e-3,fcsubexact*1e-2,\
    fcsubexact*1e-1,fcsubexact*1e0,fcsubexact*1e1,fcsubexact*1e2))
    fcexact = np.append(fcexact,1e6)

    newfc = []
    newfcexact = []
    for i in range(len(fc)):
        if fc[i] >= flims[0] and fc[i] <= flims[1]:
            newfc.append(fc[i])
            newfcexact.append(fcexact[i])
    fc = np.array(newfc)
    fcexact = np.array(newfcexact)

    step = 24/width

    fc = fc[np.arange(0,len(fc),step)]
    fcexact = fcexact[np.arange(0,len(fcexact),step)]

    df = f[1] - f[0]
    spec = []
    for i in range(len(fcexact)):
        b = 2.*width
        f1 = fcexact[i]/2.**(1./b)
        f2 = fcexact[i]*2.**(1./b)
        Qr = fcexact[i]/(f2-f1)
        Qd = (np.pi/b)/(np.sin(np.pi/b))*Qr
        Hsq = np.abs(1/(1+Qd**b*((f/fcexact[i])-(fcexact[i]/f))**b))
        spec.append(np.sum(Gxx*Hsq)*df)


    spec = np.array(spec)

    return [fc,spec]