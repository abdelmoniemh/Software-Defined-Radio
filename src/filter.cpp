/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include <cmath>
#include "iostream"

void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
	// allocate memory for the impulse response
	h.clear(); h.resize(num_taps, 1.0);

	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab

	float norm_cutoff=Fc/(Fs/2.0);


	float neum=0.0;
	float denom=0.0;

	for(int i = 0; i < num_taps; i++){
		if(i == (num_taps-1)/2) {
			h[i]=norm_cutoff;
		}
		else{
			neum=norm_cutoff*std::sin(PI*norm_cutoff*(i-(num_taps-1)/2.0));
			denom=(PI*norm_cutoff*(i-(num_taps-1.0)/2.0));
			h[i]=neum/denom;
		}
		
		h[i]=h[i]*std::sin(i*PI/num_taps)*std::sin(i*PI/num_taps);

	}
}

void BPF(float fb, float fe,float fs, unsigned short int num_taps, std::vector<float> &h)
{
	h.clear(); h.resize(num_taps, 0.0);
	float normCenter = ((fe+fb)/2.0)/(fs/2.0);
	float normPass = (fe-fb)/(fs/2.0);
	float neum,denom=0.0;


	for(int i=0;i<num_taps;i++){
		if(i==((num_taps-1.0)/2.0)){
			h[i]=normPass;
		}
		else{
			neum=normPass * (sin(PI*(normPass/2.0)*(i-(num_taps-1.0)/2.0)));
			denom=PI*(normPass/2.0)*(i-(num_taps-1.0)/2.0);

			h[i]=neum/denom;

		}

		h[i]*=cos(i*PI*normCenter);
		h[i]*=sin((i*PI)/num_taps)*sin((i*PI)/num_taps);
	}

}

// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"; this is based on
// your Python code from the take-home exercise from the previous lab
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{
	//SINGLE PASS CONVOLUTION
	y.clear(); y.resize(x.size()+h.size()-1, 0.0);

	for(int n=0;n<y.size();n++){
		for(int k=0;k<h.size();k++){
			if((n-k)>=0 && (n-k)<x.size()){
				y[n]+= h[k] * x[n-k];
			}
		}
	}

}

void convolveAndResample(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, const int up, const int down)
{
	y.clear(); y.resize((int)((up*x.size()+h.size())/down) , 0.0); 

	for(int n=0;n<y.size();n++){
		int phase = (n*down)%up;
		for(int k=phase;k<h.size();k+=up){
			int j = ((n*down-k)/up);
			if(j>=0 && j<x.size()){
				y[n]+= h[k] * x[j];
			}
			//k += up;
		}
	}

}

void blockConvolution(const std::vector<float> &x, const std::vector<float> &IR, const int &m, std::vector<float> &state, std::vector<float> &y){
	y.clear(); y.resize(x.size(), 0.0);
	for (int n = 0; n<y.size(); n++){
		y[n] = 0;
		for (int k = 0; k<IR.size(); k++){
			if (n-k>=0){
				y[n] += IR[k]*x[n-k];
			} else {
				if (m == 0) {
					y[n] += IR[k]*x[0];
				} else {
					y[n] += IR[k]*state[state.size()+n-k]; // n-k negative not indexing properly
				}
			}
		}
	}
	
	state = std::vector<float>(x.begin() + x.size() - IR.size(), x.end());
}

void convolveAndResampleBlock(const std::vector<float> &x, const std::vector<float> &IR, const int &m, std::vector<float> &state, std::vector<float> &y,
	const int up, const int down){
	y.clear(); y.resize(x.size()*((float)up/(float)down), 0.0);

	for (int n = 0; n<y.size(); n++){
		int phase = (n*down)%up;
		y[n] = 0.0;
		for (int k = phase; k<IR.size(); k+=up){
			int j = ((n*down-k)/up);
			if (j>=0){
				y[n] += IR[k]*x[j]*up;
			} else {
				if (m == 0) {
					y[n] += IR[k]*x[0]*up;
				} else {
					y[n] += IR[k]*state[state.size()+j]*up; // n-k negative not indexing properly
				}
			}
		}
	}

	state = std::vector<float>(x.begin() + x.size() - IR.size(), x.end());
	
}

void multiConvolveAndResampleBlock(
	const std::vector<float> &x,
	const std::vector<float> &z_in,
	const std::vector<float> &IR,
	const std::vector<float> &IR2,
	const int &m, 
	std::vector<float> &stateY,
	std::vector<float> &stateZ,
	std::vector<float> &y,
	std::vector<float> &z,
	const int up,
	const int down){

	y.resize(x.size()*((float)up/(float)down));
	z.resize(x.size()*((float)up/(float)down));


	for (int n = 0; n<y.size(); n++){
		int phase = (n*down)%up;
		y[n] = 0.0;
		z[n] = 0.0;
		for (int k = phase; k<IR.size(); k+=up){
			int j = ((n*down-k)/up);
			if (j>=0){
				y[n] += IR[k]*x[j]*up;
				z[n] += IR2[k]*z_in[j]*up;
			} else {
				if (m == 0) {
					y[n] += IR[k]*x[0]*up;
					z[n] += IR2[k]*z_in[0]*up;
				} else {
					y[n] += IR[k]*stateY[stateY.size()+j]*up; 
					z[n] += IR2[k]*stateZ[stateZ.size()+j]*up; 
				}
			}
		}
	}

	stateY = std::vector<float>(x.begin() + x.size() - IR.size(), x.end());
	stateZ = std::vector<float>(z_in.begin() + z_in.size() - IR.size(), z_in.end());

	
}




std::vector<float> my_lfilter(std::vector<float> &yb, const std::vector<float> &xb, const std::vector<float> &h, int block_size, std::vector<float> state){
	// allocate memory for the output (filtered) data
	yb.clear(); 
	yb.resize(xb.size(), 0.0);

	std::vector<float> b(block_size, 0.0);

	for(int n = 0; n<yb.size(); n++){	//for n in output block
		for(int k = 0; k<h.size();k++){	//for k in filter vector
			if(n-k >= 0 ){
				yb[n]+=h[k]*xb[n-k];
			}
			else{
				yb[n] = yb[n] + h[k]*state[n-k+block_size];
			}
		}
	}
	return xb;
}



void fmDemod(std::vector<float>& fm_demod, const std::vector<float>& inPhase, float& inPhaseM1, const std::vector<float>& quadPhase, float& quadPhaseM1){
	fm_demod = std::vector<float>(inPhase.size(), 0.0);
	for (int i = 0; i < inPhase.size(); i++){
		fm_demod[i] = (inPhase[i]*(quadPhase[i]-quadPhaseM1) - quadPhase[i]*(inPhase[i]-inPhaseM1))/(pow(inPhase[i],2) + pow(quadPhase[i],2));
		inPhaseM1 = inPhase[i]; quadPhaseM1 = quadPhase[i];
	}
}





std::vector<float> impulseResponseRootRaisedCosine(int Fs, int N_taps) {
    //duration for each symbol - should NOT be changed for RDS!
    double T_symbol = 1/2375.0;

    //roll-off factor (must be greater than 0 and smaller than 1)
    //for RDS a value in the range of 0.9 is a good trade-off between
    //the excess bandwidth and the size/duration of ripples in the time-domain
    double beta = 0.90;

    //the RRC impulse response that will be computed in this function
    std::vector<float> impulseResponseRRC(N_taps);

    for (int k = 0; k < N_taps; k++) {
        double t = double((k - N_taps/2.0)/ Fs);
        if (t == 0.0) {
            impulseResponseRRC[k] = double(1.0 + beta*((4 / M_PI) - 1.0));
        } 
        else if (t == -T_symbol / (4.0 * beta) || t == T_symbol / (4.0 * beta)) {
            impulseResponseRRC[k] = double((beta / sqrt(2.0)) * (((1.0 + 2.0 / M_PI) * (sin(M_PI / (4.0 * beta)))) + ((1.0 - 2.0 / M_PI) * (cos(M_PI / (4.0 * beta))))));
        } 
        else {
            impulseResponseRRC[k] = double((sin(M_PI * t * (1.0 - beta) / T_symbol) + 4.0 * beta * (t / T_symbol) * cos(M_PI * t * (1.0 + beta) / T_symbol)) / (M_PI * t * (1.0 - (4.0 * beta * t / T_symbol) * (4.0 * beta * t / T_symbol)) / T_symbol));
        }
    }

    //returns the RRC impulse response to be used by convolution
    return impulseResponseRRC;
}