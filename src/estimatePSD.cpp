#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"


std::vector<float> arange(int start, int end, float increment){
	std::vector<float> output;
	for (float i = start; i<end; i+=increment){
		output.push_back(i);
	}

	return output;
}



void estimatePSD(const std::vector<float> &samples, int Fs, std::string out_fname){

	std::vector<float> psd_list;


	int freq_bins = NFFT;	

	float df = (float)Fs/(float)freq_bins;
	float _df = 0;

	std::vector<float> freq=std::vector<float>((Fs/2)/df,0.0);


	//perform the arangement
	for(int i=0;i<(Fs/2)/df;i++){
		freq[i]=_df;
		_df=_df+df;
	}

	std::vector<float> hann(freq_bins,0.0);

	for(int i=0;i<freq_bins;i++){
		hann[i]=pow(sin((float)i*PI/(float)freq_bins),2);
	}


	int no_segments = int(floor(samples.size()/float(freq_bins)));
	
	for(int k=0;k<no_segments;k++){
		
		std::vector<float> windowed_samples(freq_bins,0.0);

		windowed_samples = std::vector<float>(samples.begin()+(k*freq_bins), samples.begin()+((k+1)*freq_bins));

		for(int i=0;i<freq_bins;i++){
			windowed_samples[i]=windowed_samples[i]*hann[i];
		}

		std::vector<std::complex<float>> Xf(freq_bins,0.0);

		DFT(windowed_samples, Xf);

		//center freq results

		Xf=std::vector<std::complex<float>>(Xf.begin(),Xf.begin()+(int(freq_bins/2)));

		std::vector<float> psd_seg(Xf.size(),0.0);

		for(int i=0;i<=Xf.size();i++){
			psd_seg[i]=(1.0/((float)Fs*(float)freq_bins/2.0)) * (float)(pow((float)abs(Xf[i]),2.0));
		}

		//psd_list.insert(std::end(psd_list), std::begin(psd_seg), std::end(psd_seg));
		for(int i=0;i<psd_seg.size();i++){
			psd_seg[i]*=2;
			psd_seg[i]=10*log10(psd_seg[i]);
			psd_list.push_back(psd_seg[i]);
		}

	}


	std::vector<float> psd_est=std::vector<float>(freq_bins/2,0.0);

	for(int k=0;k<freq_bins/2;k++){
		for(int l=0;l<no_segments;l++){
			psd_est[k]+=psd_list[k+l*int(freq_bins/2)];
		}
		psd_est[k] = psd_est[k] / no_segments;
	}


	for(int i=0;i<psd_est.size();i++){
		std::cout<<i<<": "<<psd_est[i]<<"\n";
	}

	logVector(out_fname, freq ,psd_est);

}



// void estimatePSD(const std::vector<float>& samples, const int Fs, std::string out_fname){

// 	std::vector<float> psd_est;
// 	std::vector<float> freq;
// 	const int freq_bins = 512;
// 	float df = (float)Fs/freq_bins; //check how it rounds if error

// 	freq = arange(0, Fs/2, df);

// 	std::vector<float> hann;
// 	hann.resize(freq_bins);
// 	for (int i = 0; i<hann.size(); i++){
// 		hann[i] = pow(sin(i*PI/(float)freq_bins),2);
// 	}

// 	int no_segments = floor((float)samples.size()/(float)freq_bins);

// 	std::vector<float> psd_list;

// 	for (int k = 0; k<no_segments; ++k){

// 		std::vector<float> windowed_samples(samples.begin()+k*freq_bins, samples.begin()+(k+1)*freq_bins);
// 		for (int i = 0; i<windowed_samples.size(); i++){
// 			windowed_samples[i] *= hann[i];
// 		}
// 		/*
// 		std::vector<std::complex<float>> Xf(windowed_samples.size()*2);
// 		std::vector<std::complex<float>> twiddles;
// 		compute_twiddles(twiddles);
// 		FFT_optimized(windowed_samples,Xf, twiddles);
// 		*/
// 		std::vector<std::complex<float>> Xf;
// 		DFT(windowed_samples,Xf);
// 		Xf = std::vector<std::complex<float>>(Xf.begin(), Xf.begin() + (int)freq_bins/2);
// 		for (int i = 0; i<Xf.size(); i++){
// 			float signalPower = 2 * (1/((float)Fs*(float)freq_bins/2.0f)) * pow(std::abs(Xf[i]),2);
// 			signalPower = 10*log10(signalPower);
// 			psd_list.push_back(signalPower);
// 		}
// 	}

// 	psd_est = std::vector<float>(freq_bins/2, 0.0);

// 	for (int k = 0; k<psd_est.size(); ++k){ 
// 		for (int l = 0; l<no_segments; ++l){
// 			psd_est[k] += psd_list[k + l*(int)freq_bins/2];
// 		}
// 		psd_est[k] /= no_segments;
// 	}

// 	std::ofstream outfile;
// 	outfile.open(out_fname);
// 	for (int i = 0; i<freq.size(); i++){
// 		outfile << freq[i] << ", " << psd_est[i] << std::endl;
// 	}
// 	outfile.close();
// }
