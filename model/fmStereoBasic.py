
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import wavfile
from scipy import signal, fft
from math import sin,pi,ceil
import math

import scipy


from fmStereoBlock import myBandpass,plotAndSavePSD
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


def firwinImpl(Ntaps, Fc, Fs):
	h=np.zeros(Ntaps)
	norm_cutoff= 2* (Fc/Fs)
	for i in range(Ntaps):
		if (i==(Ntaps-1)/2):
			h[i]=(norm_cutoff)
		else:
			h[i]=(norm_cutoff * (sin(pi*norm_cutoff*(i-(Ntaps-1)/2)))/(pi*norm_cutoff*(i-(Ntaps-1)/2)))

		h[i]=h[i]* (sin((i*pi)/Ntaps) * sin((i*pi)/Ntaps))

	return h

def singlePassImpl(audioData,audioCoeff):
	audioOut=np.zeros(len(audioCoeff) + len(audioData) -1)
	
	for n in range (len(audioOut)):
		for k in range (len(audioCoeff)):
			if n-k>=0 and n-k < len(audioData):
				audioOut[n]+= audioCoeff[k]*audioData[n-k]
	return audioOut


def myBandPassFirwin(freqLow, freqHigh, Fs, tapAmount):

	firwin_coeff = signal.firwin(
		tapAmount, [freqLow/(Fs/2), freqHigh/(Fs/2)], pass_zero=False)
	# freqzPlot(firwin_coeff, Fs, 'firwin for ' + str(int(Fc)) + ' Hz cutoff with ' + str(N_taps) + ' taps')
	return firwin_coeff



def plotAndSavePSD(signal,name):

	# plot PSD of selected block after FM demodulation
	ax0.clear()
	fmPlotPSD(ax0, signal, (rf_Fs/rf_decim)/1e3, subfig_height[0],
				'Demodulated FM (block ' + ')'+ " - Python ")
	# output binary file name (where samples are written from Python)
	fm_demod_fname = "data/fm_demod" + ".bin"
	# create binary file where each sample is a 32-bit float
	fm_demod = np.asarray(signal, np.float32)
	fm_demod.astype('float32').tofile(fm_demod_fname)

	# save figure to file
	fig.savefig("data/fmMonoBlock" +name+ ".png")

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
	w = w * Fs/(2*np.pi) # needed to draw the frequency on the X axis

	# plots the magnitude response where the x axis is normalized in rad/sample
	fig, ax1 = plt.subplots()
	ax1.set_title('Digital filter frequency response (' + msg + ')')
	ax1.plot(w, 20 * np.log10(abs(h)), 'b')
	ax1.set_ylabel('Amplitude [dB]', color='b')
	ax1.set_xlabel('Frequency [Hz]')



# flag that keeps track if your code is running for
# in-lab (il_vs_th = 0) vs takehome (il_vs_th = 1)

if __name__ == "__main__":

	# read the raw IQ data from the recorded file
	# IQ data is assumed to be in 8-bits unsigned (and interleaved)
	in_fname = "data/AudioIQsamples/samples5.raw"

	#in_fname = "data/AudioIQsamples/samples0.raw"


	raw_data = np.fromfile(in_fname, dtype='uint8')
	print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
	# IQ data is normalized between -1 and +1 in 32-bit float format
	iq_data = (np.float32(raw_data) - 128.0)/128.0
	print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

	# coefficients for the front-end low-pass filter
	rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))

	# filter to extract the FM channel (I samples are even, Q samples are odd)
	i_filt = signal.lfilter(rf_coeff, 1.0, iq_data[0::2])
	q_filt = signal.lfilter(rf_coeff, 1.0, iq_data[1::2])

	# downsample the FM channel
	i_ds = i_filt[::rf_decim]
	q_ds = q_filt[::rf_decim]

	# set up the subfigures for plotting
	subfig_height = np.array([1]) # relative heights of the subfigures
	plt.rc('figure', figsize=(10, 5))	# the size of the entire figure
	fig, ax0 = plt.subplots(nrows=1, gridspec_kw={'height_ratios': subfig_height})
	fig.subplots_adjust(hspace = .6)

	# FM demodulator (check the library)
	fm_demod, dummy = fmDemodArctan(i_ds, q_ds)
	# we use a dummy because there is no state for this single-pass model

	# IF DONE

	Fs = 240000.0           # sampling rate
	Fc = 16000.0            # cutoff frequency
	N_taps = 151          # number of taps for the FIR
	# coefficients for the filter to extract mono audio

	audio_coeff = signal.firwin(N_taps, Fc/(Fs/2), window=('hann'))

	# Stereo coeffs
	stereo_coeff = myBandpass(22e3, 54e3, Fs, N_taps)


	# Pilot coeffs
	pilot_coeffs = myBandpass(18.5e3, 19.5e3, Fs, N_taps)


	audio_filt = signal.lfilter(audio_coeff, 1.0, fm_demod)
	stereo_filt= signal.lfilter(stereo_coeff, 1.0, fm_demod)
	pilot_filt = signal.lfilter(pilot_coeffs, 1.0, fm_demod)

	logArray(pilot_coeffs, "pilot_coeffs")

	logArray(pilot_filt, "pilot_filt")

	pllOut,_=fmPll(pilot_filt,19e3,Fs,ncoScale=2.0)

	logArray(pllOut, "pllOut")


	delay_amount=int((N_taps-1)/2)
	delayed_mono=[0.0]*delay_amount


#	for i in range (len(audio_filt)):
#		delayed_mono.append(float(audio_filt[i]))

	#delayed_mono=(delayed_mono+audio_filt)
	delayed_mono=np.append(delayed_mono,audio_filt)

	# mixing...
	mixed_stereo=pllOut[:-1]*stereo_filt*2.0


	mixed_stereo_filtered = signal.lfilter(audio_coeff, 1.0, mixed_stereo)
	



	plotAndSavePSD(mixed_stereo,"_mixed_stereo")
	plotAndSavePSD(mixed_stereo_filtered,"_mixed_stereo_filtered")


	fb_pilot = 18.5e3
	fe_pilot = 19.5e3
	fs = 240e3
	N_taps = 151




	delayed_mono = delayed_mono[::5]
	mixed_stereo_filtered=mixed_stereo_filtered[::5]

	left_channel = (delayed_mono[:len(mixed_stereo_filtered)]-mixed_stereo_filtered)
	right_channel = (delayed_mono[:len(mixed_stereo_filtered)]+mixed_stereo_filtered)




	freqzPlot(pilot_coeffs, 240000,"my BPF impl")
	plotAndSavePSD(left_channel,"_left_channel")
	plotAndSavePSD(right_channel,"_right_channel")
	plotAndSavePSD(fm_demod,"_fm_demod")
	plotAndSavePSD(pilot_filt,"_pilot_filt")
	plotAndSavePSD(stereo_filt/2,"_stereo_filt")
	plotAndSavePSD(pllOut,"_pllOut")


	left_channel=16384.0*left_channel
	right_channel=16384.0*right_channel
	delayed_mono=16384.0*delayed_mono


	stereo_data = np.vstack((left_channel, right_channel))
	stereo_data = stereo_data.transpose()
	wavfile.write("data/StereoBasic_L_R.wav",int(audio_Fs),np.int16((stereo_data)))


	# Write mono data
	out_fname = "data/StereoBasic.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16(delayed_mono))

	# Write left channel
	out_fname = "data/StereoBasic_L.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16(left_channel))
	
	# Write right channel
	out_fname = "data/StereoBasic_R.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16(right_channel)) # (L+R)-(L-R) = 2R



	# Write stereo channel (L-R) - No mono component
	out_fname = "data/StereoBasic_Stereo_only.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16((mixed_stereo_filtered)*16384.0))

