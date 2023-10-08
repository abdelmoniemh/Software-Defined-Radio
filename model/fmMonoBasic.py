#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

"""
The command-line instructions for recording RF data are only for those
who have the RF dongle (Nooelec NESDR Smart v4 Bundle) available at home.
After you have installed the drivers to work with the RF dongle,
the 8-bit unsigned values for the I/Q pairs can be recorded as follows:

rtl_sdr -f 99.9M -s 2.4M - > iq_samples.raw

The above assumes that we are tuned to the FM station at 99.9 MHz,
we use an RF sample rate of 2.4 Msamples/sec and our file is called
iq_samples.raw (change as you see fit).

For the above use case, the data acquisition runs indefinitely,
hence the recording needs to be stopped by pressing Ctrl+C.
If we wish to stop it after a pre-defined number of samples,
e.g., 12 million I/Q pairs (5 seconds at 2.4 Msamples/sec),
we can use an extra argument:

rtl_sdr -f 99.9M -s 2.4M -n 12000000 - > iq_samples.raw

To check if the raw I/Q data has been recorded properly, place the file
in the "data" sub-folder from your project repo and run this Python file
from the "model" sub-folder. It should produce both the .png image files
(of the PSD estimates) for a couple of blocks, as well as the .wav file.

In the source code below (check lines 90+) you can observe where the
raw_data is read and the normalization of the 8-bit unsigned I/Q samples
to 32-bit float samples (in the range -1 to +1) is done; while the
32-bit floats and the range -1 to +1 are optional choices (used by
many third-party SDR software implementations), it is at the discretion
of each project group to decide how to handle the 8-bit unsigned I/Q samples
in their Python model and C++ implementation.
"""

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
# should be the same as rf_Fs / rf_decim / audio_decim


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


<<<<<<< HEAD:model/fmMonoBasic.py
=======
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

	# uncomment the lines below if you wish to inspect the phase response
	# Note: as important as the phase response is for some applications,
	# it is not critical at this stage because we expect a linear phase in the passband

	# ax2 = ax1.twinx()
	# angles = np.unwrap(np.angle(h))
	# ax2.plot(w, angles, 'g')
	# ax2.set_ylabel('Angle (radians)', color='g')

	#plt.show()


>>>>>>> stereo channel math done in cpp. Graphing scripts for cpp-> python PSD plots:model/fmStereoBasic.py
# flag that keeps track if your code is running for
# in-lab (il_vs_th = 0) vs takehome (il_vs_th = 1)
il_vs_th = 1

if __name__ == "__main__":

	# read the raw IQ data from the recorded file
	# IQ data is assumed to be in 8-bits unsigned (and interleaved)
<<<<<<< HEAD:model/fmMonoBasic.py
	in_fname = "data/AudioIQsamples/samples9.raw"
=======
	in_fname = "data/AudioIQsamples/samples1.raw"
>>>>>>> stash:model/fmStereoBasic.py

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

	# PSD after FM demodulation
	fmPlotPSD(ax0, fm_demod, (rf_Fs/rf_decim)/1e3, subfig_height[0],\
			'Demodulated FM (full recording)')

	Fs = 240000           # sampling rate
	Fc = 16000            # cutoff frequency
	N_taps = 51          # number of taps for the FIR
	# coefficients for the filter to extract mono audio
	if il_vs_th == 0:
		# to be updated by you during the in-lab session based on firwin
		# same principle  as for rf_coeff (but different arguments, of course)

<<<<<<< HEAD:model/fmMonoBasic.py
=======
	audio_coeff = signal.firwin(N_taps, Fc/(Fs/2), window=('hann'))

	# Stereo coeffs
	stereo_coeff = myBandpass(22e3, 54e3, Fs, N_taps)

	# Pilot coeffs
	pilot_coeffs = myBandpass(18.5e3, 19.5e3, Fs, N_taps)

	audio_filt = signal.lfilter(audio_coeff, 1.0, fm_demod)

	stereo_filt= signal.lfilter(stereo_coeff, 1.0, fm_demod)
	pilot_filt = signal.lfilter(pilot_coeffs, 1.0, fm_demod)


	logArray(pilot_filt, "pilot_filt")


	pllOut=fmPll(pilot_filt,19e3,Fs,2.0)

	logArray(pllOut, "pllOut")


	pllOut=np.array(pllOut)
	stereo_filt=np.array(stereo_filt)

	mixed_stereo=pllOut[:-1]*stereo_filt

	#mixed_stereo=signal.lfilter(pllOut,1.0,stereo_filt)*2
>>>>>>> stereo channel math done in cpp. Graphing scripts for cpp-> python PSD plots:model/fmStereoBasic.py

		#firwin_coeff = signal.firwin(N_taps, Fc/(Fs/2), window=('hann'))


		audio_coeff = signal.firwin(N_taps, Fc/(Fs/2), window=('hann'))
	else:
		# to be updated by you for the takehome exercise
		# with your own code for impulse response generation
		audio_coeff = firwinImpl(N_taps, Fc, Fs)

<<<<<<< HEAD:model/fmMonoBasic.py
	# extract the mono audio data through filtering
	if il_vs_th == 0:
		# to be updated by you during the in-lab session based on lfilter
		# same principle as for i_filt or q_filt (but different arguments)
		#filteredI = signal.lfilter(firwin_coeff, 1.0, i_ds)
		
		audio_filt = signal.lfilter(audio_coeff, 1.0, fm_demod)
	else:
		# to be updated by you for the takehome exercise
		# with your own code for single pass convolution
		audio_filt = singlePassImpl(fm_demod, audio_coeff)
=======


	plotAndSavePSD(mixed_stereo,"_mixed_stereo")
	plotAndSavePSD(mixed_stereo_filtered,"_mixed_stereo_filtered")
>>>>>>> stash:model/fmStereoBasic.py


	fb_pilot = 18.5e3
	fe_pilot = 19.5e3
	fs = 240e3
	N_taps = 151

	# you should uncomment the plots below once you have processed the data

	# PSD after extracting mono audio
	# fmPlotPSD(ax1, audio_filt, (rf_Fs/rf_decim)/1e3, subfig_height[1], 'Extracted Mono')

	# downsample audio data (see the principle for i_ds or q_ds)
	audio_data = audio_filt[::5] # to be updated by you during in-lab (same code for takehome)

	# PSD after decimating mono audio
	# fmPlotPSD(ax2, audio_data, audio_Fs/1e3, subfig_height[2], 'Downsampled Mono Audio')

<<<<<<< HEAD:model/fmMonoBasic.py
<<<<<<< HEAD:model/fmMonoBasic.py
	# save PSD plots
	#fig.savefig("../data/fmMonoBasic.png")
	#plt.show()
=======
	left_channel = (audio_data/2+mixed_stereo_filtered/2)
	right_channel = (audio_data/2-mixed_stereo_filtered/2)

	freqzPlot(pilot_coeffs, 240000,"my BPF impl")

=======
	# devide by 2.0 after adding/subtracting
	left_channel = (mixed_stereo_filtered+audio_data)/2.0
	right_channel = (mixed_stereo_filtered-audio_data)/2.0
>>>>>>> stash:model/fmStereoBasic.py


	logArray(audio_data, "audio_data")
	logArray(mixed_stereo_filtered, "mixed_stereo_filtered")
	freqzPlot(pilot_coeffs, 240000,"my BPF impl")
	plotAndSavePSD(left_channel,"_left_channel")
	plotAndSavePSD(right_channel,"_right_channel")
	plotAndSavePSD(fm_demod,"_fm_demod")
	plotAndSavePSD(pilot_filt,"_pilot_filt")
	plotAndSavePSD(stereo_filt/2,"_stereo_filt")
	plotAndSavePSD(pllOut,"_pllOut")


	#left_channel=left_channel
	#right_channel=right_channel

<<<<<<< HEAD:model/fmMonoBasic.py
	stereo_data = np.vstack((left_channel, right_channel))
	stereo_data = stereo_data.transpose()
	wavfile.write("data/test.wav",int(audio_Fs),np.int16((stereo_data/2)))
>>>>>>> stereo channel math done in cpp. Graphing scripts for cpp-> python PSD plots:model/fmStereoBasic.py

	# write audio data to file (assumes audio_data samples are -1 to +1)
	out_fname = "data/fmMonoBasic.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16((audio_data/2)*32767))
=======

	# write audio data to file (assumes audio_data samples are -1 to +1)

	# check possible issue with scaling
	out_fname = "data/StereoBasic.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16((audio_data)*16384.0))	# Mono Block ONLY


	out_fname = "data/StereoBasic_L.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16(left_channel*16384.0))	# (L+R)+(L-R) = 2L


	out_fname = "data/StereoBasic_R.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16(right_channel*16384.0)) # (L+R)-(L-R) = 2R


	#Stereo_audio_final=np.concatenate((left_channel*16384.0, right_channel*16384.0), axis=0)


	out_fname = "data/StereoBasic_Stereo_only.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16(mixed_stereo_filtered*16384.0))





>>>>>>> stash:model/fmStereoBasic.py
	# during FM transmission audio samples in the mono channel will contain
	# the sum of the left and right audio channels; hence, we first
	# divide by two the audio sample value and then we rescale to fit
	# in the range offered by 16-bit signed int representation
	print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")
