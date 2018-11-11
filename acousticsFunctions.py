"""
acousticsFunctions includes binfileload and weighting functions


"""

import numpy as np
import struct
import sys

def binfileload(path, IDname, IDnum, CHnum, N=10, NStart=0):
    IDnum = "%03.0f" %IDnum
    CHnum = "%03.0f" %CHnum
    if sys.platform.startswith('win'):
        filename = path+"\\"+IDname+IDnum+"_"+CHnum+".bin"
    else:
        filename = path+"/"+IDname+IDnum+"_"+CHnum+".bin"
    N = int(N)
    with open(filename,'rb') as fin:
        num_data_bytes = 4*N
        data_str = fin.read(num_data_bytes)
        fmt = str(N)+'f'
        data = struct.unpack(fmt,data_str)  #N 4-byte floats
    return np.array(data)





"""
(W,Gain) = weighting(f,type='A')
This function returns the weighting curves, W evaluated at the frequencies,
f. Valid types are 'A','B','C','D','E','G','U','ITUR468', and 'M'.  If type is not specified, the
default is A-weighting.  To apply the weighting function to a power or
autospectrum, the spectrum is multiplied by this function, W.. 
Sources: Wikipedia (A-weighting) and https://en.wikipedia.org/wiki/ITU-R_468_noise_weighting
Author: Kent Gee 
"""
from numpy import sqrt

def weighting(f,type='A'):
    type = type.upper()
    if type == 'A':
        K = 10.**(2./20.)
        W = K*(12200.**2*f**4)/(f**2+20.6**2)/(f**2+12200.**2)/sqrt(f**2+107.7**2)/sqrt(f**2+737.9**2)
        W = W**2
    elif type == 'B':
        K=10**(.17/20)
        W=K*(12200**2*f**3)/(f**2+20.6**2)/(f**2+12200.**2)/sqrt(f**2+158.5**2)
        W=W**2
    elif type == 'C':
        K=10**(.06/20)
        W=K*(12200.**2*f**2)/(f**2+20.6**2)/(f**2+12200.**2)
        W=W**2
    elif type == 'Ds':
        pass
    elif type == 'E':
        pass
    else:
        print('Unknown weighting type... more to come')
    
    return (W, 10*np.log10(W))


'''
switch type
    case 'A' %wikipedia.  Matches ANSI standard
        K=10^(2/20);
        W=K*(12200^2*f.^4)./(f.^2+20.6^2)./(f.^2+12200^2)./sqrt(f.^2+107.7^2)./sqrt(f.^2+737.9^2);
        W=W.^2;
    case 'B' %wikipedia.  Matches ANSI standard
        K=10^(.17/20);
        W=K*(12200^2*f.^3)./(f.^2+20.6^2)./(f.^2+12200^2)./sqrt(f.^2+158.5^2);
        W=W.^2;
    case 'C' %wikipedia. Matches ANSI standard
        K=10^(.06/20);
        W=K*(12200^2*f.^2)./(f.^2+20.6^2)./(f.^2+12200^2);
        W=W.^2;
    case 'Ds' s-domain expression from Wikipedia  
        K=91104.32;  -domain expression
        s=j*2*pi*f;
        W=abs(K*s.*(s.^2+6532*s+4.0975e7)./(s+1776.3)./(s+7288.5)./(s.^2+21514*s+3.8836e8));
        W=W.^2;
    case 'D' f-domain expression from ANSI standard
        K=2.1024164e8;
        W=K*f.^2.*((-519.8)^2+(f+876.2).^2).*((-519.8)^2+(f-876.2).^2)./((-282)^2+f.^2)./((-1160)^2+f.^2)./((-1712)^2+(f+2628).^2)./((-1712)^2+(f-2628).^2);
    case 'E' % ANSI standard
        K=3.8341500e16;
        W=K*f.^4.*((-735)^2+(f+918).^2).*((-735)^2+(f-918).^2)./((-53.5)^2+f.^2)./((-378)^2+f.^2)./((-865)^2+f.^2)./((-4024)^2+(f+3966).^2)./((-4024)^2+(f-3966).^2)./((-6500)^2+f.^2);
    case 'G' %ANSI standard
        K=10^(231.992/20); determined empirically. Matches table B2 but not Fig. 31, which has the wrong gain
        W=K*f.^8./((-.707)^2+(f+.707).^2)./((-.707)^2+(f-.707).^2)./((-19.27)^2+(f+5.16).^2)./((-19.27)^2+(f-5.16).^2)./((-14.11)^2+(f+14.11).^2)./((-14.11)^2+(f-14.11).^2)./((-5.16)^2+(f+19.27).^2)./((-5.16)^2+(f-19.27).^2);
    case 'U' %ANSI standar
        K=10^(490.183/10); etermined empirically
        W=K./((-12200)^2+(f).^2).^2./((-7850)^2+(f+8800).^2)./((-7850)^2+(f-8800).^2)./((-2900)^2+(f+12150).^2)./((-2900)^2+(f-12150).^2);
    case 'ITUR468' %wikipedia 
        K=10^(18.2/20);
        h1=-4.737338981378384e-24*f.^6+2.043828333606125e-15*f.^4-1.363894795463638e-7*f.^2+1;
        h2=1.306612257412824e-19*f.^5-2.118150887518656e-11*f.^3+5.559488023498642e-4*f;
        W=K*1.246332637532143e-4*f./sqrt(h1.^2+h2.^2);
        W=W.^2;
    case 'M' %wikipedia 468 article referencing ISO standard.  ALso referred to as CCIR/ARM 
        K=10^(18.2/20)*10^(-5.5905/20); %puts 0 dB at 2 kHz instead of 1 kHz.
        h1=-4.737338981378384e-24*f.^6+2.043828333606125e-15*f.^4-1.363894795463638e-7*f.^2+1;
        h2=1.306612257412824e-19*f.^5-2.118150887518656e-11*f.^3+5.559488023498642e-4*f;
        W=K*1.246332637532143e-4*f./sqrt(h1.^2+h2.^2);
        W=W.^2;
    otherwise
        disp('Unknown weighting type');
end
Gain=10*log10(W);
end
'''