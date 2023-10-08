#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import numpy as np
import math

import matplotlib.pyplot as plt



# integrator = 0.0
# phaseEst = 0.0
# feedbackI = 1.0
# feedbackQ = 0.0
# ncoOut[0] = 1.0
# trigOffset = 0	



def fmPll(pllIn, freq, Fs, ncoScale = 1.0, phaseAdjust = 0.0, normBandwidth = 0.01, state=[0.0,0.0,1.0,0.0,1.0,0.0]):


	# scale factors for proportional/integrator terms
	# these scale factors were derived assuming the following:
	# damping factor of 0.707 (1 over square root of 2)
	# there is no oscillator gain and no phase detector gain
	Cp = 2.666
	Ci = 3.555

	# gain for the proportional term
	Kp = (normBandwidth)*Cp
	# gain for the integrator term
	Ki = (normBandwidth*normBandwidth)*Ci

	# output array for the NCO
	ncoOut = np.empty(len(pllIn)+1)

	# initialize internal state
	integrator = state[0]
	phaseEst = state[1]
	feedbackI = state[2]
	feedbackQ = state[3]
	ncoOut[0] = state[4]
	trigOffset = state[5]
	# note: state saving will be needed for block processing

	errorI_log=[]
	errorQ_log=[]

	for k in range(len(pllIn)):

		# phase detector
		errorI = pllIn[k] * (+feedbackI)  # complex conjugate of the
		errorQ = pllIn[k] * (-feedbackQ)  # feedback complex exponential

		errorI_log.append(errorI)
		errorQ_log.append(errorQ)


		# four-quadrant arctangent discriminator for phase error detection
		errorD = math.atan2(errorQ, errorI)

		# loop filter
		integrator = integrator + Ki*errorD

		# update phase estimate
		phaseEst = phaseEst + Kp*errorD + integrator

		# internal oscillator
		trigOffset += 1
		trigArg = 2*math.pi*(freq/Fs)*(trigOffset) + phaseEst
		feedbackI = math.cos(trigArg)
		feedbackQ = math.sin(trigArg)
		ncoOut[k+1] = math.cos(trigArg*ncoScale + phaseAdjust)

	# for stereo only the in-phase NCO component should be returned
	# for block processing you should also return the state

		state[0] = integrator
		state[1] = phaseEst
		state[2] = feedbackI
		state[3] = feedbackQ
		state[4] = ncoOut[0]
		state[5] = trigOffset


	return ncoOut, state
	# for RDS add also the quadrature NCO component to the output

if __name__ == "__main__":

	pass
