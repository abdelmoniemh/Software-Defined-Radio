import matplotlib.pyplot as plt
import numpy as np
from scipy.io import wavfile
from scipy import signal, fft
from math import sin, pi, ceil
import math

import scipy


#from fmStereoBlock import myBandpass, plotAndSavePSD
from PLL.fmPll import fmPll

# use fmDemodArctan and fmPlotPSD
from fmSupportLib import fmDemodArctan, fmPlotPSD
# for take-home add your functions

# the radio-frequency (RF) sampling rate
# this sampling rate is either configured on RF hardware
# or documented when a raw file with IQ samples is provided
rf_Fs = 2.4e6

# the cutoff frequency to extract the FM channel from raw IQ data
rf_Fc = 100e3

# the number of taps for the low-pass filter to extract the FM channel
# this default value for the width of the impulse response should be changed
# depending on some target objectives, like the width of the transition band
# and/or the minimum expected attenuation from the pass to the stop band
rf_taps = 151

# the decimation rate when reducing the front end sampling rate (i.e., RF)
# to a smaller samping rate at the intermediate frequency (IF) where
# the demodulated data will be split into the mono/stereo/radio data channels
rf_decim = 10

# audio sampling rate (we assume audio will be at 48 KSamples/sec)
audio_Fs = 48e3
# should be the same ∫as rf_Fs / rf_decim / audio_decim





P=[
    [1,0,0,0,0,0,0,0,0,0],           #0
    [0,1,0,0,0,0,0,0,0,0],           #1
    [0,0,1,0,0,0,0,0,0,0],           #2
    [0,0,0,1,0,0,0,0,0,0],           #3
    [0,0,0,0,1,0,0,0,0,0],           #4
    [0,0,0,0,0,1,0,0,0,0],           #5
    [0,0,0,0,0,0,1,0,0,0],           #6
    [0,0,0,0,0,0,0,1,0,0],           #7
    [0,0,0,0,0,0,0,0,1,0],           #8
    [0,0,0,0,0,0,0,0,0,1],           #9
    [1,0,1,1,0,1,1,1,0,0],           #10
    [0,1,0,1,1,0,1,1,1,0],           #11
    [0,0,1,0,1,1,0,1,1,1],           #12
    [1,0,1,0,0,0,0,1,1,1],           #13
    [1,1,1,0,0,1,1,1,1,1],           #14
    [1,1,0,0,0,1,0,0,1,1],           #15
    [1,1,0,1,0,1,0,1,0,1],           #16
    [1,1,0,1,1,1,0,1,1,0],           #17
    [0,1,1,0,1,1,1,0,1,1],           #18
    [1,0,0,0,0,0,0,0,0,1],           #19
    [1,1,1,1,0,1,1,1,0,0],           #20
    [0,1,1,1,1,0,1,1,1,0],           #21
    [0,0,1,1,1,1,0,1,1,1],           #22
    [1,0,1,0,1,0,0,1,1,1],           #23
    [1,1,1,0,0,0,1,1,1,1],           #24
    [1,1,0,0,0,1,1,0,1,1]            #25
]



Syndromes={
      "A" :[1,1,1,1,0,1,1,0,0,0],
      "B" :[1,1,1,1,0,1,0,1,0,0],
      "C" :[1,0,0,1,0,1,1,1,0,0],
      "C'":[1,1,1,1,0,0,1,1,0,0],
      "D" :[1,0,0,1,0,1,1,0,0,0]
}




def firwinImpl(Ntaps, Fc, Fs):
	h = np.zeros(Ntaps)
	norm_cutoff = 2 * (Fc/Fs)
	for i in range(Ntaps):
		if (i == (Ntaps-1)/2):
			h[i] = (norm_cutoff)
		else:
			h[i] = (norm_cutoff * (sin(pi*norm_cutoff*(i-(Ntaps-1)/2))) /
			        (pi*norm_cutoff*(i-(Ntaps-1)/2)))

		h[i] = h[i] * (sin((i*pi)/Ntaps) * sin((i*pi)/Ntaps))

	return h


def singlePassImpl(audioData, audioCoeff):
	audioOut = np.zeros(len(audioCoeff) + len(audioData) - 1)

	for n in range(len(audioOut)):
		for k in range(len(audioCoeff)):
			if n-k >= 0 and n-k < len(audioData):
				audioOut[n] += audioCoeff[k]*audioData[n-k]
	return audioOut


def myBandPassFirwin(freqLow, freqHigh, Fs, tapAmount):

	firwin_coeff = signal.firwin(
		tapAmount, [freqLow/(Fs/2), freqHigh/(Fs/2)], pass_zero=False)
	# freqzPlot(firwin_coeff, Fs, 'firwin for ' + str(int(Fc)) + ' Hz cutoff with ' + str(N_taps) + ' taps')
	return firwin_coeff


def plotAndSavePSD(signal, name):

	# plot PSD of selected block after FM demodulation
	ax0.clear()
	fmPlotPSD(ax0, signal, (rf_Fs/rf_decim)/1e3, subfig_height[0],
				'Demodulated FM (block ' + ')' + " - Python ")
	# output binary file name (where samples are written from Python)
	fm_demod_fname = "data/fm_demod" + ".bin"
	# create binary file where each sample is a 32-bit float
	fm_demod = np.asarray(signal, np.float32)
	fm_demod.astype('float32').tofile(fm_demod_fname)

	# save figure to file
	fig.savefig("data/fmMonoBlock" + name + ".png")

	return


def logArray(array, arrayName):
	fileName = "logs/"+arrayName+".txt"
	with open(fileName, "w") as file:
		for i in range(len(array)):
			line = "["+str(i)+"]"+": "+str(array[i])
			file.write(f"{line}\n")


def myBandpass(fb, fe, fs, Ntaps):

	normCenter = ((fe+fb)/2.0)/(fs/2.0)
	normPass = (fe-fb)/(fs/2.0)

	h = [0.0]*Ntaps

	for i in range(Ntaps):
		if (i == ((Ntaps-1.0)/2.0)):
			h[i] = normPass

		else:
			neum = normPass * \
				math.sin(math.pi*(normPass/2.0)*(i-(Ntaps-1.0)/2.0))
			denom = math.pi*(normPass/2.0)*(i-(Ntaps-1)/2.0)
			h[i] = neum/denom
		h[i] = h[i]*math.cos(i*math.pi*normCenter)
		h[i] = h[i]*(math.sin(i*math.pi/Ntaps)**2)

	return h


def freqzPlot(coeff, Fs, msg):

	# find the frequency response using freqz from SciPy:
	# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.freqz.html
	w, h = signal.freqz(coeff)

	# Reminder: np.pi rad/sample is actually the Nyquist frequency
	w = w * Fs/(2*np.pi)  # needed to draw the frequency on the X axis

	# plots the magnitude response where the x axis is normalized in rad/sample
	fig, ax1 = plt.subplots()
	ax1.set_title('Digital filter frequency response (' + msg + ')')
	ax1.plot(w, 20 * np.log10(abs(h)), 'b')
	ax1.set_ylabel('Amplitude [dB]', color='b')
	ax1.set_xlabel('Frequency [Hz]')


def conv_resampler_fast(x,h,m,state,du,ds): #du is U ds is N
    #y=np.zeros(int((du*len(x)+len(h))//ds))
    y=np.zeros(int((len(x)*(du/ds))))

    for n in range(len(y)):
        #define phase for each value of y we want to calculate
        phase=int(((n*ds)%du))
        #k=phase
        for k in range(phase,len(h)):
            
             
            #concolve x by h
            j= int( ((n*ds-k)/du) )

            #convolve
            if(j>=0):
                y[n]+=h[k]*x[j] #Must use J for x here
            else:
                if(m==0):
                    y[n] += h[k]*x[0]
                else:
                    y[n] += h[k]*state[len(state)+j]
                    
                    
            k+=du

            
           # k+=du #must increment k by this much
    state=x[-(len(h)-1):]
    return y,state

def conv_downsample_fast(x,h,ds): #The fast way to compute convolution with downsampling
    y=np.zeros((1,len(x)+len(h))/2)
    for n in range(y):
        for k in range(h):
            N=ds*N
            if n-k>=0 and n-k <len(x):
                y[n]+= h[k]*x[N]
    return y

def conv_resampler_slow(x,h,du,ds):
    ex = np.zeros(du*len(x)) #longer version of x
    for i in range(len(x)):#this loop is for zero padding
        ex[i*du]=x[i] #spacing out our input values so we have add the zero padding

    y=conv_downsample_fast(x,h,ds) #this will do the convolution and downsampling together

    return y

def fmDemodImpl(inPhase, inPhaseM1, quadPhase, quadPhaseM1):
	m = np.zeros(len(inPhase))
	for i in range(len(inPhase)):
		m[i] = (inPhase[i]*(quadPhase[i]-quadPhaseM1) - quadPhase[i]*(inPhase[i]-inPhaseM1))/(inPhase[i]**2 + quadPhase[i]**2)
		inPhaseM1, quadPhaseM1 = inPhase[i], quadPhase[i]
	return m, inPhaseM1, quadPhaseM1

def impulseResponseRootRaisedCosine(Fs, N_taps):

	"""
	Root raised cosine (RRC) filter
	Fs  		sampling rate at the output of the resampler in the RDS path
				sampling rate must be an integer multipler of 2375
				this integer multiple is the number of samples per symbol
	N_taps  	number of filter taps
	"""

	# duration for each symbol - should NOT be changed for RDS!
	T_symbol = 1/2375.0

	# roll-off factor (must be greater than 0 and smaller than 1)
	# for RDS a value in the range of 0.9 is a good trade-off between
	# the excess bandwidth and the size/duration of ripples in the time-domain
	beta = 0.90

	# the RRC inpulse response that will be computed in this function
	impulseResponseRRC = np.empty(N_taps)

	for k in range(N_taps):
		t = float((k-N_taps/2))/Fs
		# we ignore the 1/T_symbol scale factor
		if t == 0.0: impulseResponseRRC[k] = 1.0 + beta*((4/math.pi)-1)
		elif t == -T_symbol/(4*beta) or t == T_symbol/(4*beta):
			impulseResponseRRC[k] = (beta/np.sqrt(2))*(((1+2/math.pi)* \
					(math.sin(math.pi/(4*beta)))) + ((1-2/math.pi)*(math.cos(math.pi/(4*beta)))))
		else: impulseResponseRRC[k] = (math.sin(math.pi*t*(1-beta)/T_symbol) +  \
					4*beta*(t/T_symbol)*math.cos(math.pi*t*(1+beta)/T_symbol))/ \
					(math.pi*t*(1-(4*beta*t/T_symbol)*(4*beta*t/T_symbol))/T_symbol)

	# returns the RRC impulse response to be used by convolution
	return impulseResponseRRC

def findLocalMaxMin(data, symbol_rate):
    max = -1000000
    max_i = 0
    min = 1000000
    min_i = 0
    for i in range(symbol_rate):
        if data[i] > max:
            max_i = i 
            max = data[i]
        if data[i] < min:
            min_i = i
            min = data[i]
    if abs(max) > abs(min):
        return max_i
    else:
        return min_i 



def matchSyndrome(block):
    blocks=["A","B","C'","C", "D"]
    for OffsetType in blocks:
        if(block==Syndromes[OffsetType]):
            return OffsetType
    return -1

def getSyndrome(block):
    matrix_col=[0]*26
    result=[0]*10

    for i in range(10):
        for j in range(26):
            matrix_col[j]=P[j][i]
            
        prod = np.array(block,dtype=bool) & np.array(matrix_col,dtype=bool)
        
        for elem in prod:
            result[i]^=elem

        #result[i] = '0b1' if prod.count(1) % 2 else '0b0'

    return result


if __name__ == "__main__":


    in_fname = "data/AudioIQsamples/iq_samples.raw"

    raw_data = np.fromfile(in_fname, dtype='uint8')
    print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
    # IQ data is normalized between -1 and +1 in 32-bit float format
    iq_data = (np.float32(raw_data) - 128.0)/128.0
    print("Reformatted raw RF data to 32-bit float format (" +
        str(iq_data.size * iq_data.itemsize) + " bytes)")

    # coefficients for the front-end low-pass filter
    rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))

    # filter to extract the FM channel (I samples are even, Q samples are odd)
    i_filt = signal.lfilter(rf_coeff, 1.0, iq_data[0::2])
    q_filt = signal.lfilter(rf_coeff, 1.0, iq_data[1::2])
    print("I/Q Filtered")


    # downsample the FM channel
    i_ds = i_filt[::rf_decim]
    q_ds = q_filt[::rf_decim]

    fm_demod, dummy = fmDemodArctan(i_ds, q_ds)
    print("FM Demodulated")


###################### RF FRONT END DONE #########################


##################### STARTING RDS PROCESSING ########################



    Fs = 240000.0           # sampling rate
    Fc = 16000.0            # cutoff frequency
    N_taps = 151          # number of taps for the FIR


    # FILTER RDS SIGNAL 
    rds_coeff = myBandPassFirwin(54e3, 60e3, Fs, N_taps)
    rds_filt = signal.lfilter(rds_coeff, 1.0, fm_demod)
    print("Filtered RDS")
    rds_symbol_rate = 11



    # SQUARE RDS SIGNAL & FILTER TO GENERATE CARRIER
    rds_carrier = rds_filt*rds_filt
    rds_carrier_coeff = myBandPassFirwin(113.5e3, 114.5e3, Fs, N_taps)
    rds_carrier_filt = signal.lfilter(rds_carrier_coeff, 1.0, rds_carrier)
    print("Filtered RDS Carrier")


    pllOut, _ =fmPll(rds_carrier_filt,114e3,Fs,ncoScale=0.5) # Double check ncoScale
    rds_carrier_i = pllOut[:-1]
    print("Pll Done")


    # PASS FILTERED RDS THROUGH ALLPASS BEFORE MIXING FOR DELAY
    all_pass_coeff = np.zeros(N_taps)
    all_pass_coeff[(N_taps-1)//2] = 1
    rds_delayed = signal.lfilter(all_pass_coeff, 1.0, rds_filt)
    print("Delayed rds_filt")


    # MIXING
    rds_mixed_i = 2*rds_delayed*rds_carrier_i
#    rds_mixed_q = 2*rds_delayed*rds_carrier_q   # FOR DEBUGGING ONLY





    # FILTER MIXED VALUE AT 3KHZ
    rds_demod_coeff = signal.firwin(N_taps, 3e3/(Fs/2), window=('hann')) #Low Pass filter 3kHz
    rds_mixed_filt_i = signal.lfilter(rds_demod_coeff, 1.0, rds_mixed_i)
#    rds_mixed_filt_q = signal.lfilter(rds_demod_coeff, 1.0, rds_mixed_q)


    Up = 209
    Down = 1920

    # Relational Resampler
    rds_resampled_i = signal.resample_poly(rds_mixed_filt_i, Up, Down)
    rds_resampled_i *=Up


#    rds_resampled_q = signal.resample_poly(rds_mixed_filt_q, Up, Down)
#    rds_resampled_q *=Up

    
    # RRC FILTERING 
    cosine_coeff = impulseResponseRootRaisedCosine(11*2375, N_taps)
    rds_demod_i = signal.lfilter(cosine_coeff, 1.0, rds_resampled_i)
#    rds_demod_q = signal.lfilter(cosine_coeff, 1.0, rds_resampled_q)
    
    plt.plot(range(440),rds_demod_i[:440])
    plt.savefig("rds_cosine_i.png")
    plt.clf()


    delay = int((N_taps)*Up/Down) + (N_taps-1)//2 
    rds_demod_i = rds_demod_i[delay:]
#    rds_demod_q = rds_demod_q[delay:]

 

    start_i = findLocalMaxMin(rds_demod_i, 11)
    start_q = findLocalMaxMin(rds_demod_q, 11)


    #rds_demod_i = rds_demod_i[start_i::11]
    #rds_demod_q = rds_demod_q[start_q::11]
    print(len(rds_demod_i), len(rds_demod_q))
    plt.scatter( rds_demod_i[:1000]/max(rds_demod_i[:1000]) ,rds_demod_q[:1000])
    plt.savefig("rds_demod.png")
    plt.clf()

    plt.plot(range(440),rds_demod_i[:440])
    plt.savefig("rds_demod_i.png")
    plt.clf()
    plt.plot(range(440),rds_demod_q[:440])
    plt.savefig("rds_demod_q.png")
    plt.clf()
    

def decoding(toBeDecoded):
    i=0
    cdr = []
    while(i<len(toBeDecoded)-1):
        if((toBeDecoded[i]>0 and toBeDecoded[i+1]>0) or (toBeDecoded[i]<0 and toBeDecoded[i+1]<0)):
            #cdr=[]
            #print(f"clear cdr at i = {i} because of {toBeDecoded[i], toBeDecoded[i+1]}")
            #i+=1
            #toBeDecoded=toBeDecoded[i+1:]
            #print(f"remaining input {toBeDecoded}")
            #i=0
            cdr.append(0)
            i+=2
            continue
        if (i+1 >= len(toBeDecoded)):
              break 
        if(toBeDecoded[i]<toBeDecoded[i+1]):
            cdr.append(0)
            i+=2
        else:
            cdr.append(1)
            i+=2

    post_cdr=[]
    #print(f"before xor = {cdr}")
    for x in range(len(cdr)-1):
        post_cdr.append(cdr[x] ^ cdr[x+1])

    return post_cdr

#print(f"demod to 100 {rds_demod_i[:100]}")

block=[]
cdrOut = decoding(rds_demod_i)


for i in cdrOut:
      block.append(i)
      if(len(block)==26):
            syndrome=getSyndrome(block)
            offsetType=matchSyndrome(syndrome)

            if(offsetType=='A'):
                a = int("".join(str(i) for i in block[:15]),2)
                print("A:", hex(a))
                #print(block[:15])
            elif(offsetType=='B'):
                print("B: ",block[:15])
            elif(offsetType=='C'):
                print("C: ",block[:15])
            elif(offsetType=='C\''):
                print("C': ",block[:15])
            elif(offsetType=='D'):
                print("D: ",block[:15])
                                          

            if(offsetType!=-1):
                #print(offsetType)
                #print(block[:15])
                block=[]
            else:
                block.pop(0)


#print(f"output: {cdrOut}")

#print(f"output: {cdrOut}")
