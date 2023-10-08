#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
from math import sin,pi,ceil
import math

# use fmDemodArctan and fmPlotPSD
from fmSupportLib import fmDemodArctan, fmPlotPSD
# for take-home add your functions

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
	h=np.zeros(Ntaps)
	norm_cutoff= 2* (Fc/Fs)
	for i in range(Ntaps):
		if (i==(Ntaps-1)/2):
			h[i]=(norm_cutoff)
		else:
			h[i]=(norm_cutoff * (sin(pi*norm_cutoff*(i-(Ntaps-1)/2)))/(pi*norm_cutoff*(i-(Ntaps-1)/2)))

		h[i]=h[i]* (sin((i*pi)/Ntaps) * sin((i*pi)/Ntaps))

	return h

def blockFilteringImpl(audioData,audioCoeff,state,m):
	audioOut=np.zeros(len(audioData))
	for n in range (len(audioOut)):
		for k in range (len(audioCoeff)):
			if n-k>=0:
				audioOut[n]+= audioCoeff[k]*audioData[n-k]
			else:
				if m==0:
					audioOut[n]+= audioCoeff[k]*audioData[0]
				else:
					audioOut[n]+= audioCoeff[k]*state[n-k]
	state=audioData[-(len(audioCoeff)-1):]
				
	return audioOut,state

def fmDemodImpl(inPhase, inPhaseM1, quadPhase, quadPhaseM1):
	m = np.zeros(len(inPhase))
	for i in range(len(inPhase)):
		m[i] = (inPhase[i]*(quadPhase[i]-quadPhaseM1) - quadPhase[i]*(inPhase[i]-inPhaseM1))/(inPhase[i]**2 + quadPhase[i]**2)
		inPhaseM1, quadPhaseM1 = inPhase[i], quadPhase[i]
	return m, inPhaseM1, quadPhaseM1


def logArray(array,arrayName):
	fileName="logs/"+arrayName+".txt"
	with open(fileName, "w") as file:
		for i in range(len(array)):
			line="["+str(i)+"]"+": "+str(array[i])
			file.write(f"{line}\n")


def conv_resampler_fast(x,h,du,ds): #du is U ds is N
    y=np.zeros((du*len(x)+len(h))/ds)

    for n in range(len(y)):
        #define phase for each value of y we want to calculate
        phase=((n*ds)%du)
        #k=phase
        for k in range(phase,len(h)):
            j= int(((n*ds-k)/du))
            if(j>=0 and j<len(x)):
                y[n]+=h[k]*x[j] #Must use J for x here

    return y

if __name__ == "__main__":

	# read the raw IQ data from the recorded file
	# IQ data is assumed to be in 8-bits unsigned (and interleaved)
	in_fname = "data/AudioIQsamples/iq_samples.raw"
	raw_data = np.fromfile(in_fname, dtype='uint8')
	print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
	# IQ data is normalized between -1 and +1 in 32-bit float format
	iq_data = (np.float32(raw_data) - 128.0)/128.0
	print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

	# coefficients for the front-end low-pass filter
	rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))

	Fs = 240000          # sampling rate
	Fc = 16000            # cutoff frequency
	N_taps = 151          # number of taps for the FIR
	# coefficients for the filter to extract mono audio
	if il_vs_th == 0:
		# to be updated by you during the in-lab session based on firwin
		# same principle  as for rf_coeff (but different arguments, of course)

		audio_coeff = signal.firwin(N_taps, Fc/(Fs/2), window=('hann'))
	else:
		# to be updated by you for the takehome exercise
		# with your own code for impulse response generation
		audio_coeff = firwinImpl(N_taps, Fc, Fs)

	# set up the subfigures for plotting
	subfig_height = np.array([0.8, 2, 1.6]) # relative heights of the subfigures
	plt.rc('figure', figsize=(7.5, 7.5))	# the size of the entire figure
	fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
	fig.subplots_adjust(hspace = .6)

	# select a block_size that is a multiple of KB
	# and a multiple of decimation factors
	block_size = 1024 * rf_decim * audio_decim * 2
	block_count = 0

	# states needed for continuity in block processing
	state_i_lpf_100k = np.zeros(rf_taps-1)
	state_q_lpf_100k = np.zeros(rf_taps-1)
	state_phase = 0

	
	# add state as needed for the mono channel filter
	filter_state = np.zeros(N_taps-1)

	# audio buffer that stores all the audio blocks
	audio_data = np.array([]) # used to concatenate filtered blocks (audio data)

	inPhaseM1 = 0
	quadPhaseM1 = 0

	# if the number of samples in the last block is less than the block size
	# it is fine to ignore the last few samples from the raw IQ file
	while (block_count+1)*block_size < len(iq_data):

		# if you wish to have shorter runtimes while troubleshooting
		# you can control the above loop exit condition as you see fit
		print('Processing block ' + str(block_count))

		# filter to extract the FM channel (I samples are even, Q samples are odd)
		i_filt, state_i_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
				iq_data[(block_count)*block_size:(block_count+1)*block_size:2],
				zi=state_i_lpf_100k)
		q_filt, state_q_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
				iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2],
				zi=state_q_lpf_100k)
		# downsample the I/Q data from the FM channel
		i_ds = i_filt[::rf_decim]
		q_ds = q_filt[::rf_decim]

		# FM demodulator
		if il_vs_th == 0:
			# already given to you for the in-lab
			# take particular notice of the "special" state-saving
			fm_demod, state_phase = fmDemodArctan(i_ds, q_ds, state_phase)
		else:
			# you will need to implement your own FM demodulation based on:
			# https://www.embedded.com/dsp-tricks-frequency-demodulation-algorithms/
			# see more comments on fmSupportLib.py - take particular notice that
			# you MUST have also "custom" state-saving for your own FM demodulator
			fm_demod, inPhaseM1, quadPhaseM1 = fmDemodImpl(i_ds, inPhaseM1, q_ds, quadPhaseM1)

		# extract the mono audio data through filtering
		if il_vs_th == 0:
			audio_filt, filter_state = signal.lfilter(audio_coeff, 1.0, fm_demod, zi=filter_state)
		else:
			audio_filt, filter_state = blockFilteringImpl(fm_demod, audio_coeff, filter_state, block_count)
		# 	# to be updated by you for the takehome exercise
		# 	# with your own code for BLOCK convolution
		#print(len(audio_filt), len(filter_state), len(fm_demod), len(audio_coeff))


		# downsample audio datas
		# to be updated by you during in-lab (same code for takehome)
		audio_block = audio_filt[::5]

		# concatenate the most recently processed audio_block
		# to the previous blocks stored already in audio_data
		#
		audio_data = np.concatenate((audio_data, audio_block))
		#

		# to save runtime select the range of blocks to log data
		# this includes both saving binary files as well plotting PSD
		# below we assume we want to plot for graphs for blocks 10 and 11
		if block_count >= 10 and block_count < 12:

			# plot PSD of selected block after FM demodulation
			ax0.clear()
			fmPlotPSD(ax0, fm_demod, (rf_Fs/rf_decim)/1e3, subfig_height[0], \
					'Demodulated FM (block ' + str(block_count) + ')')
			# output binary file name (where samples are written from Python)
			fm_demod_fname = "data/fm_demod_" + str(block_count) + ".bin"
			# create binary file where each sample is a 32-bit float
			fm_demod.astype('float32').tofile(fm_demod_fname)

			# plot PSD of selected block after extracting mono audio
			# ... change as needed

			# plot PSD of selected block after downsampling mono audio
			# ... change as needed

			# save figure to file
			fig.savefig("data/fmMonoBlock" + str(block_count) + ".png")

		if(block_count==0):
			# log data for block 
			logArray(audio_coeff,"audio_coeff")
			logArray(rf_coeff,"rf_coeff")
			logArray(i_filt,"i_filt")
			logArray(q_filt,"q_filt")
			logArray(iq_data[(block_count)*block_size:(block_count+1)*block_size:2],"i_data")

			logArray(fm_demod,"fm_demod")
			logArray(audio_block,"audio_block")




		block_count += 1


	print('Finished processing all the blocks from the recorded I/Q samples')

	# write audio data to file
	print(f"len of audio data = {len(audio_data)}")
	out_fname = "data/fmMonoBlock.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16((audio_data/2)*32767))
	logArray((audio_data/2)*32767,"audio_data/2")

	print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")

	# uncomment assuming you wish to show some plots
	# plt.show()
