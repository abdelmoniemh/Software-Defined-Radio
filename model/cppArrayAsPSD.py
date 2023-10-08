
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import wavfile
from scipy import signal
from math import sin,pi,ceil
import math

import scipy

from fmStereoBlock import myBandpass,plotAndSavePSD
from PLL.fmPll import fmPll

# use fmDemodArctan and fmPlotPSD
from fmSupportLib import fmDemodArctan, fmPlotPSD



def readTextfileAsArray(filename_in):

    result=[];

    file1 = open(filename_in, 'r')
    count = 0

    while True:
        count += 1
    
        # Get next line from file
        line = file1.readline()
    
        # if line is empty
        # end of file is reached

        for i in range(len(line)):
            if(line[i]==" "):
                a=(line[i+1:-1])
                result.append(float(a))

        if not line:
            break
   
         

    return np.array(result,dtype=float)

rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 151
rf_decim = 10

audio_Fs = 48e3
audio_decim = 5
block_count=0;

def plotAndSavePSD(signal,name):

	# set up the subfigures for plotting
	subfig_height = np.array([1]) # relative heights of the subfigures
	plt.rc('figure', figsize=(10, 5))	# the size of the entire figure
	fig, (ax0) = plt.subplots(nrows=1, gridspec_kw={'height_ratios': subfig_height})
	fig.subplots_adjust(hspace = .6)


	# plot PSD of selected block after FM demodulation
	ax0.clear()
	fmPlotPSD(ax0, signal, (rf_Fs/rf_decim)/1e3, subfig_height[0],
				'Demodulated FM (block ' + str(block_count) + ') '+ " - CPP")


	# output binary file name (where samples are written from Python)

	# create binary file where each sample is a 32-bit float

	# save figure to file
	fig.savefig("../data/CPP_fmMonoBlock" + str(block_count) +name+ ".png")

	return

def plotTextfileAsPSD(filename_in, filename_out):
    plotAndSavePSD(readTextfileAsArray(filename_in),filename_out)
    return

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

	fig.savefig("../data/CPP_fmMonoBlock" +msg+"- CPP"+ ".png")


if __name__=="__main__":

    # set up the subfigures for plotting
    subfig_height = np.array([1]) # relative heights of the subfigures
    plt.rc('figure', figsize=(10, 5))	# the size of the entire figure
    fig, ax0 = plt.subplots(nrows=1, gridspec_kw={'height_ratios': subfig_height})
    fig.subplots_adjust(hspace = .6)

    freqzPlot(readTextfileAsArray("../cpp_logs/_pilot_coeffs.dat"), 240000, "cpp BPF (pilot)");    
    freqzPlot(readTextfileAsArray("../cpp_logs/_stereo_coeffs.dat"), 240000, "cpp BPF (stereo)");    


    plotTextfileAsPSD("../cpp_logs/_pilotData.dat","_pilotData")
    plotTextfileAsPSD("../cpp_logs/_fm_demod.dat","_fm_demod")
    #plotTextfileAsPSD("../cpp_logs/_stereoData.dat","_stereoData")
    plotTextfileAsPSD("../cpp_logs/_pllOut.dat","_pllOut")
    plotTextfileAsPSD("../cpp_logs/_mixed_signal.dat","_mixed_signal")
    plotTextfileAsPSD("../cpp_logs/_mixed_signal_filtered.dat","_mixed_signal_filtered")
    plotTextfileAsPSD("../cpp_logs/_monoData.dat","_monoData")




