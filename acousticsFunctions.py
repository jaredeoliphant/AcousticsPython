# modules acousticsFunctions.py includes binfileload, weighting and beamform functions
import numpy as np
import struct
import sys

def binfileload(path, IDname, IDnum, CHnum, N=10, NStart=0):
    """
    "binfileload" is used to input binary data from a file specified at a certain path with an
    ID number and an Channel number
    N number of data points needs to be specified currently. (Default is only 10 data points)
    NStart does not currently have any functionality currently either
    translated to python by Jared Oliphant
    """
    
    # format the IDnum and CHnum strings
    IDnum = "%03.0f" %IDnum
    CHnum = "%03.0f" %CHnum

    # check if the system is windows or linux
    if sys.platform.startswith('win'):
        filename = path+"\\"+IDname+IDnum+"_"+CHnum+".bin"
    else:
        filename = path+"/"+IDname+IDnum+"_"+CHnum+".bin"

    # coerce to an integer
    N = int(N)
    print('opening ',filename)
    # use struct.unpack method to get read the data string
    with open(filename,'rb') as fin:
        num_data_bytes = 4*N
        data_str = fin.read(num_data_bytes)
        fmt = '<'+str(N)+'f'
        # little-endian (@ and = give the same result)
        data = struct.unpack(fmt,data_str)  #N 4-byte floats
    
    #f = open(filename,'rb')
    #data1 = np.fromfile(f,dtype=float,count=N,sep="")
    #print(data1)
    # return as an array
    return np.array(data)















def weighting(f,type='A'):
    """
    W,Gain = weighting(f,type='A')
    Gain = 10*log10(W)
    This function returns the weighting curves, W evaluated at the frequencies,
    f. Valid types are 'A','B','C','D','E','G','U','ITUR468', and 'M'.  If type is not specified, the
    default is A-weighting.  To apply the weighting function to a power or
    autospectrum, the spectrum is multiplied by this function, W.. 
    Sources: Wikipedia (A-weighting) and https://en.wikipedia.org/wiki/ITU-R_468_noise_weighting
    Author: Kent Gee    
    translated to python by Jared Oliphant
    """

    # coerce type to upper case
    type = type.upper()

    # calculate based on the type of weighting desired
    if type == 'A':
        K = 10.0**(2/20.0)
        W = K*(12200.**2*f**4)/(f**2+20.6**2)/(f**2+12200.**2)/np.sqrt(f**2+107.7**2)/np.sqrt(f**2+737.9**2)
        W = W**2
    elif type == 'B':
        K=10**(.17/20)
        W=K*(12200**2*f**3)/(f**2+20.6**2)/(f**2+12200.**2)/np.sqrt(f**2+158.5**2)
        W=W**2
    elif type == 'C':
        K=10**(.06/20)
        W=K*(12200.**2*f**2)/(f**2+20.6**2)/(f**2+12200.**2)
        W=W**2
    elif type == 'Ds':
        K=91104.32
        s=1j*2*np.pi*f
        W=np.abs(K*s*(s**2+6532.*s+4.0975e7)/(s+1776.3)/(s+7288.5)/(s**2+21514.*s+3.8836e8))
        W=W**2
    elif type == 'D':
        K=2.1024164e8
        W=K*f**2*((-519.8)**2+(f+876.2)**2)*((-519.8)**2+(f-876.2)**2)/((-282.)**2+f**2)/((-1160.)**2+f**2)/((-1712.)**2+(f+2628.)**2)/((-1712)**2+(f-2628.)**2)
    elif type == 'E':
        K=3.8341500e16
        W=K*f**4*((-735.)**2+(f+918.)**2)*((-735.)**2+(f-918.)**2)/((-53.5)**2+f**2)/((-378.)**2+f**2)/((-865.)**2+f**2)/((-4024.)**2+(f+3966.)**2)/((-4024.)**2+(f-3966.)**2)/((-6500.)**2+f**2)
    elif type == 'G':
        K=10**(231.992/20)
        W=K*f**8/((-.707)**2+(f+.707)**2)/((-.707)**2+(f-.707)**2)/((-19.27)**2+(f+5.16)**2)/((-19.27)**2+(f-5.16)**2)/((-14.11)**2+(f+14.11)**2)/((-14.11)**2+(f-14.11)**2)/((-5.16)**2+(f+19.27)**2)/((-5.16)**2+(f-19.27)**2)
    elif type == 'U':
        K=10**(490.183/10)
        W=K/((-12200)**2+(f)**2)**2/((-7850.)**2+(f+8800.)**2)/((-7850.)**2+(f-8800.)**2)/((-2900.)**2+(f+12150.)**2)/((-2900.)**2+(f-12150.)**2)
    elif type == 'ITUR468':
        K=10**(18.2/20)
        h1=-4.737338981378384e-24*f**6+2.043828333606125e-15*f**4-1.363894795463638e-7*f**2+1
        h2=1.306612257412824e-19*f**5-2.118150887518656e-11*f**3+5.559488023498642e-4*f
        W=K*1.246332637532143e-4*f/np.sqrt(h1**2+h2**2)
        W=W**2
    elif type == 'M':
        K=10**(18.2/20)*10**(-5.5905/20)
        h1=-4.737338981378384e-24*f**6+2.043828333606125e-15*f**4-1.363894795463638e-7*f**2+1
        h2=1.306612257412824e-19*f**5-2.118150887518656e-11*f**3+5.559488023498642e-4*f
        W=K*1.246332637532143e-4*f/np.sqrt(h1**2+h2**2)
        W=W**2
    else:
        print('Unknown weighting type')

    # return weighting and the gain values (dB)
    return W, 10*np.log10(W)
























def beamform(x, fs, d, lookangles, desiredangle=None):
    """
    x is a list signals   x = [x1(t),x2(t),x3(t),...]
    fs is the sampling frequency of the signals
    d is the spacing between microphones, in meters
    angle is an array of look angles, in radians
    The function assumes a speed of sound of 343 m/s
    """
    if desiredangle == None:
        waveform = None

    ## speed of sound
    c = 343.0

     # number of signals in the list x
    channels = len(x)   

    # maximum delay required (end fire incidence)
    maxTau = (channels-1)*d/c   # (channels-1)*d is the max distance
    dt = 1.0/fs   # time between samples
    maxindices = maxTau/dt  # number of indices to reach max delay

    # take the middle chunk of data so we never go too far when we delay
    indices = int(round(maxindices)+18)
    
    # never shift the first signal
    unshifted = x[0][indices:-indices]

    print('beamforming...'," "*(len(lookangles)-14),"| end")
    SPLvec = []

    for angle in lookangles:
        theta = angle - 90.0*np.pi/180.0    # shift it 90 degrees to get into my reference (90 deg is broadside)
        shift = []

        # delay
        for j in range(channels-1):
            tau = (j+1)*d*np.sin(theta)/c    # shift j+1 down
            shiftindices = int(round(tau/dt))
            shift.append(x[j+1][indices-shiftindices:-indices-shiftindices])

        # sum / mean
        m = unshifted/float(channels)
        for i in range(len(shift)):
            m += shift[i]/float(channels)

        # return the desired waveform    
        if desiredangle != None:
            if angle == desiredangle*np.pi/180:
                waveform = m

        pref = 2e-5
        SPLvec.append(20*np.log10(np.mean(m**2)/pref**2))
        print("#",end="", flush=True)


    SPLvec = np.array(SPLvec)

    return SPLvec,waveform
