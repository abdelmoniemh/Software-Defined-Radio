#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
import matplotlib
from scipy.io import wavfile
from scipy import signal
import numpy as np
from math import sin, pi, ceil
import math

import time

from PLL.fmPll import fmPll




# use fmDemodArctan and fmPlotPSD
from fmSupportLib import fmDemodArctan, fmPlotPSD
# for take-home add your functions

matplotlib.use("GTK4Agg")

rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 151
rf_decim = 10

audio_Fs = 48e3
audio_decim = 5
# add other settings for audio, like filter taps, ...

# flag that keeps track if your code is running for
# in-lab (il_vs_th = 0) vs takehome (il_vs_th = 1)
il_vs_th = 1


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


def blockFilteringImpl(audioData, audioCoeff, state, m):
	audioOut = np.zeros(len(audioData))
	for n in range(len(audioOut)):
		for k in range(len(audioCoeff)):
			if n-k >= 0:
				audioOut[n] += audioCoeff[k]*audioData[n-k]
			else:
				if m == 0:
					audioOut[n] += audioCoeff[k]*audioData[0]
				else:
					audioOut[n] += audioCoeff[k]*state[n-k]
	state = audioData[-(len(audioCoeff)-1):]

	return audioOut, state


def fmDemodImpl(inPhase, inPhaseM1, quadPhase, quadPhaseM1):
	m = np.zeros(len(inPhase))
	for i in range(len(inPhase)):
		m[i] = (inPhase[i]*(quadPhase[i]-quadPhaseM1) - quadPhase[i] *
				(inPhase[i]-inPhaseM1))/(inPhase[i]**2 + quadPhase[i]**2)
		inPhaseM1, quadPhaseM1 = inPhase[i], quadPhase[i]
	return m, inPhaseM1, quadPhaseM1


def logArray(array, arrayName):
	fileName = "logs/"+arrayName+".txt"
	with open(fileName, "w") as file:
		for i in range(len(array)):
			line = "["+str(i)+"]"+": "+str(array[i])
			file.write(f"{line}\n")


# =====================

def generateSin(Fs, interval, frequency=7.0, amplitude=5.0, phase=0.0):

	dt = 1.0/Fs                          # sampling period (increment in time)
	time = np.arange(0, interval, dt)    # time vector over interval

	# generate the sin signal
	x = amplitude*np.sin(2*math.pi*frequency*time+phase)

	return time, x


def generateMultitoneSin():

	time, wave1 = generateSin(1024, 1.0, 100, 5, 0)
	_, wave2 = generateSin(1024, 1.0, 200, 10, 0)
	_, wave3 = generateSin(1024, 1.0, 300, 15, 0)

	wave = wave1+wave2+wave3

	plotTime(wave, time)

	return time, wave


def plotSpectrum(x, Fs, type='FFT'):

	n = len(x)             # length of the signal
	df = Fs/n              # frequency increment (width of freq bin)

	# compute Fourier transform, its magnitude and normalize it before plotting
	if type == 'FFT':
		Xfreq = np.fft.fft(x)

	if type == 'DFT':
		# C# = zeros((2,2),dtype=complex)
		Xfreq = np.zeros(n, dtype=complex)
		for m in range(n):
			for k in range(n):
				Xfreq[m] = Xfreq[m]+x[k] * np.exp(-2*1j*np.pi*k*m/n)

	XMag = abs(Xfreq)/n

	# Note: because x is real, we keep only the positive half of the spectrum
	# Note also: half of the energy is in the negative half (not plotted)
	XMag = XMag[0:int(n/2)]

	# freq vector up to Nyquist freq (half of the sample rate)
	freq = np.arange(0, Fs/2, df)

	fig, ax = plt.subplots()
	ax.plot(freq, XMag)
	ax.set(xlabel='Frequency (Hz)', ylabel='Magnitude',
		   title=type + ' Frequency domain plot')
	#fig.savefig(type + ' ' + "freq.png")
	plt.show()

	return Xfreq, XMag  # return x, y pair


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

	# uncomment the lines below if you wish to inspect the phase response
	# Note: as important as the phase response is for some applications,
	# it is not critical at this stage because we expect a linear phase in the passband

	# ax2 = ax1.twinx()
	# angles = np.unwrap(np.angle(h))
	# ax2.plot(w, angles, 'g')
	# ax2.set_ylabel('Angle (radians)', color='g')


def plotTime(y, time):

	fig, ax = plt.subplots()
	ax.plot(time, y)
	ax.set(xlabel='Time (sec)', ylabel='Amplitude',
		   title='Time domain plot')
	fig.savefig("time.png")
	plt.show()


def myBandpass(fb, fe, fs, Ntaps):

	normCenter = ((fe+fb)/2.0)/(fs/2.0)
	normPass = (fe-fb)/(fs/2.0)

	h = [0.0]*Ntaps

	for i in range(Ntaps):
		if (i == ((Ntaps-1.0)/2.0)):
			h[i] = normPass
		else:
			neum = normPass * math.sin(math.pi*(normPass/2.0)*(i-(Ntaps-1.0)/2.0))
			denom = math.pi*(normPass/2.0)*(i-(Ntaps-1)/2.0)
			h[i] = neum/denom
		h[i] = h[i]*math.cos(i*math.pi*normCenter)
		h[i] = h[i]*(math.sin(i*math.pi/Ntaps)**2)

	return h


def myBandPassFirwin(freqLow, freqHigh, Fs, tapAmount):

	firwin_coeff = signal.firwin(
		tapAmount, [freqLow/(Fs/2), freqHigh/(Fs/2)], pass_zero=False)
	# freqzPlot(firwin_coeff, Fs, 'firwin for ' + str(int(Fc)) + ' Hz cutoff with ' + str(N_taps) + ' taps')
	return firwin_coeff


def plotAndSavePSD(signal,name):



	# plot PSD of selected block after FM demodulation
	ax0.clear()
	fmPlotPSD(ax0, signal, (rf_Fs/rf_decim)/1e3, subfig_height[0],
				'Demodulated FM (block ' + str(block_count) + ')')
	# output binary file name (where samples are written from Python)
	fm_demod_fname = "data/fm_demod_" + \
		str(block_count) + ".bin"
	# create binary file where each sample is a 32-bit float
	fm_demod = np.asarray(signal, np.float32)
	fm_demod.astype('float32').tofile(fm_demod_fname)

	# plot PSD of selected block after extracting mono audio
	ax1.clear()
	fmPlotPSD(ax1, audio_filt, (rf_Fs/rf_decim)/1e3,
				subfig_height[1], 'Undownsampled Mono Audio')


	# plot PSD of selected block after downsampling mono audio
	ax2.clear()
	fmPlotPSD(ax2, audio_data, audio_Fs/1e3,
				subfig_height[2], 'Downsampled Mono Audio')

	# save figure to file
	fig.savefig("data/fmMonoBlock" + str(block_count) +name+ ".png")

	return



if __name__ == "__main__":

	fb_pilot = 18.5e3
	fe_pilot = 19.5e3
	fs = 240e3
	N_taps = 151

	fb_stereo = 22e3
	fe_stereo = 54e3

	audio_data_stereo=[0.0]





	print('Take-home exercise for the digital filter design')


	in_fname = "data/AudioIQsamples/stereo_l0_r9.raw"
	#in_fname = "data/AudioIQsamples/labroom_samples102_9.raw"

	#in_fname = "data/AudioIQsamples/samples0.raw"

	raw_data = np.fromfile(in_fname, dtype='uint8')
	print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
	# IQ data is normalized between -1 and +1 in 32-bit float format
	iq_data = (np.float32(raw_data) - 128.0)/128.0
	print("Reformatted raw RF data to 32-bit float format (" +
		   str(iq_data.size * iq_data.itemsize) + " bytes)")

	



	audio_coeff = signal.firwin(N_taps, 16000/(240000/2), window=('hann'))
	carrier_coeffs = myBandPassFirwin(fb_pilot, fe_pilot, fs, N_taps)
	stereo_coeffs = myBandPassFirwin(fb_stereo, fe_stereo, fs, N_taps)





	rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))


	delay_coeffs=[0.0]*151
	delay_coeffs[75]=1.0

	monoQ=[]

	for i in range(75):
		monoQ.append(0.0)


	subfig_height = np.array([0.8, 2, 1.6])
	plt.rc('figure', figsize=(7.5, 7.5))  # the size of the entire figure
	fig, (ax0, ax1, ax2) = plt.subplots(
		nrows=3, gridspec_kw={'height_ratios': subfig_height})
	fig.subplots_adjust(hspace=.6)


	block_size = 1024 * rf_decim * audio_decim * 2
	block_count = 0
	block_number = len(iq_data) // block_size

	state_i_lpf_100k = np.zeros(rf_taps-1)
	state_q_lpf_100k = np.zeros(rf_taps-1)
	state_phase = 0


	audio_data = np.array([])
	left_channel=np.array([])
	right_channel=np.array([])
	stereo_data=np.array([])
	audio_state = [0.0] * 150
	carrier_state = [0.0] * 150
	stereo_state = [0.0] * 150
	mixed_stereo_state =[0.0] * 150




	pllState=[0.0,0.0,1.0,0.0,1.0,0.0]


	mono_delayed_state=[0.0]*150


	while (block_count+1)*block_size < len(iq_data):


		# if you wish to have shorter runtimes while troubleshooting
		# you can control the above loop exit condition as you see fit
		print('Processing block ' + str(block_count))

		# filter to extract the FM channel (I samples are even, Q samples are odd)
		i_filt, state_i_lpf_100k = signal.lfilter(rf_coeff, 1.0,
													iq_data[(
														block_count)*block_size:(block_count+1)*block_size:2],
													zi=state_i_lpf_100k)
		q_filt, state_q_lpf_100k = signal.lfilter(rf_coeff, 1.0,
													iq_data[(
														block_count)*block_size+1:(block_count+1)*block_size:2],
													zi=state_q_lpf_100k)

		# downsample the I/Q data from the FM channel
		i_ds = i_filt[::rf_decim]
		q_ds = q_filt[::rf_decim]


		fm_demod, state_phase = fmDemodArctan(i_ds, q_ds, state_phase)


		# MONO PATH
		audio_filt, audio_state = signal.lfilter(audio_coeff, 1.0, fm_demod, zi=audio_state)						#MONO EXTRACTION
		mono_delayed=[0.0]*len(audio_filt)
		#mono_delayed, mono_delayed_state = signal.lfilter(delay_coeffs, 1.0, audio_filt, zi=mono_delayed_state)	#MONO DELAY

		for i in audio_filt:
			monoQ.append(i)

		for i in range(len(audio_filt)):
			mono_delayed[i]=monoQ.pop(0)


		# STEREO PATH
		carrier,carrier_state=signal.lfilter(carrier_coeffs, 1.0, fm_demod, zi=carrier_state)					#CARRIER EXTRACTION
		stereo,stereo_state=signal.lfilter(stereo_coeffs, 1.0, fm_demod, zi=stereo_state)						#STEREO EXTRACTION
		ncoOut,pllState=fmPll(pllIn=carrier,freq=19_000,Fs=240e3,ncoScale=2.0, state=pllState) 					#PLL 
		mixed_signal=stereo*ncoOut[0:-1]*2.0	
																		#MIXING
		mixed_stereo,mixed_stereo_state=signal.lfilter(audio_coeff, 1.0, mixed_signal, zi=mixed_stereo_state)	#MIXED EXTRACTION


		




		audio_block = mono_delayed[::5]
		audio_block_stereo_mixed = mixed_stereo[::5]


		left_channel_block=(audio_block-audio_block_stereo_mixed)
		right_channel_block=(audio_block+audio_block_stereo_mixed)


		audio_data = np.concatenate((audio_data, audio_block))
		stereo_data= np.concatenate((stereo_data, audio_block_stereo_mixed))
		left_channel = np.concatenate((left_channel, left_channel_block))
		right_channel = np.concatenate((right_channel, right_channel_block))

	
		if(block_count==1):
			logArray(ncoOut,"_pllOutPY")
			logArray(mono_delayed,"_mono_delayedPY")
			logArray(carrier,"_carrierPY")



		block_count += 1

	print('Finished processing all the blocks from the recorded I/Q samples')

	# write audio data to file

	out_fname = "data/fmMonoBlock.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16((audio_data/2)*32767))

	out_fname = "data/fmMonoBlock_stereo.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16((stereo_data)*32767))


	out_fname = "data/fmMonoBlock_stereo_mono_L.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16((left_channel/2.0)*32767.0))

	out_fname = "data/fmMonoBlock_stereo_mono_R.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16((right_channel/2.0)*32767.0))


