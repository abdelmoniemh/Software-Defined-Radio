/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include <project.h>
#include <mutex>
#include <thread>
#include <queue>
#include <memory>
#include <condition_variable>
#include <atomic>
#include "global.h"

#define QUEUE_ELEMS 10

void audio_path(Parameters &parameters,std::queue<std::shared_ptr<AudioBlock>> &queue , std::mutex &my_mutex,std::condition_variable &my_cvar, std::atomic<bool> &flag){//fucntion for audio path (mono and stereo)  thread

	std::vector<float> audio_coeff;

	std::vector<float> allData;
	std::vector<float> allDataLeft;
	std::vector<float> allDataRight;

	
	float audio_fc = 16000;
	int maxNumOfTaps = (((float)parameters.blockSize/2.0)/(float)parameters.rf_down) - 1;
	unsigned int num_taps_audio = 151*parameters.if_up > maxNumOfTaps ? maxNumOfTaps : 151*parameters.if_up;
	//unsigned int num_taps_audio = maxNumOfTaps;
	std::cerr << "max num of taps " << maxNumOfTaps << std::endl;
	std::cerr << "actual num of taps " << num_taps_audio << std::endl;

	
	std::vector<float> filterState(num_taps_audio-1, 0.0);
	
	impulseResponseLPF(parameters.if_fs * parameters.if_up, audio_fc, num_taps_audio, audio_coeff);
	
	PLL pll((((float)parameters.blockSize)/(float)parameters.rf_down));

	std::vector<float> stateCarrier(num_taps_audio-1, 0.0);	//make this more optimized ...
	std::vector<float> stateStereo(num_taps_audio-1, 0.0);	//make this more optimized ...
	std::vector<float> delayState(parameters.blockSize, 0.0);	//make this more optimized ...
	std::vector<float> stateStereoMixed(num_taps_audio, 0.0);	//make this more optimized ...

	std::vector<float> delay_coeffs(num_taps_audio, 0.0);
	delay_coeffs[((num_taps_audio-1)/2.0)*((float)parameters.if_up/(float)parameters.if_down)]=1.0;

	std::vector<float> pilot_coeff;
	BPF(parameters.carrier_lower_freq, parameters.carrier_upper_freq, parameters.if_fs, 151, pilot_coeff);

	std::vector<float> stereo_coeff;
	BPF(parameters.stereo_lower_freq, parameters.stereo_upper_freq, parameters.if_fs, 151, stereo_coeff);


	std::queue<float> monoQ;	// QUEUE FOR THE MONO AUDIO DELAY

	//int delay_amount=(((151-1)/2.0))*(float)parameters.if_up/(float)parameters.if_down + 1;
	int delay_amount = 75; //optimal delay value
	std::cerr << "Delay Amount  = " << delay_amount << std::endl;
	for(int i=0;i<delay_amount;i++){
        monoQ.push(0.0);    //add some delay to mono
    }


	while (true) { //loop to process all block by poping RF proccessed blocks from queue
		std::unique_lock<std::mutex> my_lock(my_mutex);//mutex lock 
		while (queue.empty()){
        	my_cvar.wait(my_lock);
        }         
		
		std::shared_ptr<AudioBlock> block = queue.front();
		queue.pop();
		my_cvar.notify_one();
		my_lock.unlock();
		std::cerr << "Read Block " << block->getBlockId() << " from queue" << std::endl;


		block->processMonoAudio( //method to process the mono audio which will be passed onto stereo
			audio_coeff,
			delay_coeffs,
			filterState,
			delayState,
			monoQ,
			parameters.if_up,
			parameters.if_down
		);


		if (parameters.stereo) { //if stereo mode is selected, process stereo audio
			block->processStereoAudio(
				audio_coeff,
				pilot_coeff, 
				stereo_coeff,

				stateCarrier,
				stateStereo,
				stateStereoMixed,
				
				parameters.if_up,
				parameters.if_down,
				parameters.if_fs,

				pll
			);



			block->outputAudio(parameters.stereo, true, &allDataLeft, &allDataRight);
			
		} else {	
			block->outputAudio(parameters.stereo, true, &allData);
		}

		if (flag == true) //flag that is triggered in RF front end thread to indicate all blocks are read and processed
		{
			std::cerr << "Flag is One" << std::endl;
			if (queue.empty()) {//wait until queue is empty meaning all blocks have been processed in this thread
				std::cerr << "Queue is Empty" << std::endl;
				std::string out_fname = "../data/output/sample1.bin";
				if (parameters.stereo){
					//write_audio_data(out_fname, allDataLeft, allDataRight);
					write_audio_data(out_fname, allDataLeft, allDataRight);
				}
				else
					write_audio_data(out_fname, allData, allData);
				break;
			}
		}
	}
}


void rf_frontend(Parameters &parameters, std::queue<std::shared_ptr<AudioBlock>> &fm_queue,std::queue<std::shared_ptr<AudioBlock>> &fm_queue_rds, std::mutex &my_mutex,std::mutex &my_mutex_rds, std::condition_variable &my_cvar,std::condition_variable &my_cvar_rds, std::atomic<bool> &flag) //fucntion for rf front end thread
{

	unsigned int num_taps_rf = 151;

	float RF_FcUD = (parameters.if_fs / 2.0) * ((float)parameters.if_up / (float)parameters.if_down);
	float RF_Fc = RF_FcUD > parameters.if_fs / 2.0 ? parameters.if_fs / 2.0 : RF_FcUD;

	std::cerr << "Began Generating Filter Coefficiencts" << std::endl;
	std::vector<float> rf_coeff;
	impulseResponseLPF(parameters.rf_fs, 100e3, num_taps_rf, rf_coeff);
	
	std::cerr << "Done Generating Filter Coefficiencts" << std::endl;

	std::vector<float> stateI(num_taps_rf, 0.0);
	std::vector<float> stateQ(num_taps_rf, 0.0);

	std::vector<float> audio_data;
	float inPhaseM1 = 0;
	float quadPhaseM1 = 0;


	std::cerr << "Block Size of " << parameters.blockSize << std::endl;
	for (int block_id = 0;; block_id++)
	{
		std::vector<float> block_data(parameters.blockSize, 0.0);
		readStdinBlockData(parameters.blockSize, block_id, block_data);

		if ((std::cin.rdstate()) != 0)
		{
			flag = true;
			std::cerr << "End of input stream reached" << std::endl;
			break;
		}

		std::shared_ptr<AudioBlock> block = std::make_shared<AudioBlock>(block_data, block_id, parameters.blockSize);

		std::cerr << "Read block " << block_id << std::endl;

		block->preprocessRFData(
			stateI,
			stateQ,
			inPhaseM1,
			quadPhaseM1,
			rf_coeff,
			parameters.rf_up,
			parameters.rf_down);

		std::unique_lock<std::mutex> my_lock(my_mutex);

		while (fm_queue.size() >= QUEUE_ELEMS)
		{
			my_cvar.wait(my_lock);
		}
		fm_queue.push(block);
		my_cvar.notify_one();
		my_lock.unlock();

		if(parameters.rds){

		std::unique_lock<std::mutex> my_lock1(my_mutex_rds);

		//Uncomment this when we acctually use the rds thread

		

		

		while (fm_queue_rds.size() >= QUEUE_ELEMS)
		{
			my_cvar_rds.wait(my_lock1);
		}
		fm_queue_rds.push(block);
		my_cvar_rds.notify_one();
		my_lock1.unlock();

		}

	}
}

void rdbs_path(Parameters &parameters, std::queue<std::shared_ptr<AudioBlock>> &queue, std::mutex &my_mutex, std::condition_variable &my_cvar, std::atomic<bool> &flag){//fucntion for rds thread

	std::vector<float> rds_coeff;
	std::vector<float> rds_carrier_coeff;
	std::vector<float> rds_demod_coeff;
	std::vector<float> cosine_coeff;


	// FILTER DEFS.

	BPF(54e3, 60e3, parameters.if_fs, 151, rds_coeff); 					//BPF coeff for RDS extraction (actual rds extraction from radio wave)
	BPF(113.5e3, 114.5e3, parameters.if_fs, 151, rds_carrier_coeff);	//BPF coeff for RDS carrier extraction
	impulseResponseLPF(parameters.if_fs*parameters.rdsUp, 3e3, 151, rds_demod_coeff); 	//Filter for after PLL
	cosine_coeff=impulseResponseRootRaisedCosine(parameters.sps*2375, 151); 	//RRC Filter coeffs


	// STATE DEFS
	
	std::vector<float> rds_filt_state(150,0.0);			// RDS 
	std::vector<float> rds_carrier_filt_state(150,0.0);
	std::vector<float> mixed_rds_state(150,0.0);
	std::vector<float> rds_demod_i_state(150,0.0);

	std::vector<bool> manchester_state;
	std::vector<bool> diff_state;
	std::vector<int> frame_state;


	// MISC
	int start_i_state=0;


	// RDS PLL
	PLL rds_pll((((float)parameters.blockSize)/(float)parameters.rf_down));

	std::queue<float> rds_Q;

	for(int i=0;i<75;i++){
		rds_Q.push(0.0);	// ADD DELAY
	}

	while(true){
		std::unique_lock<std::mutex> my_lock(my_mutex);
		while (queue.empty()){
        	my_cvar.wait(my_lock);
        }         
		
		std::shared_ptr<AudioBlock> rdbs_block = queue.front();
		queue.pop();
		my_cvar.notify_one();
		my_lock.unlock();
		std::cerr << "Read Block " << rdbs_block->getBlockId() << " from queue RBDS" << std::endl;

		//Block processing for rds here
		rdbs_block->processRBDS(    
			//Filters
			rds_coeff,
			rds_carrier_coeff,
			rds_demod_coeff,
			cosine_coeff,

			//DELAY QUEUE,
			rds_Q,

			//Filter States
			rds_filt_state,
			rds_carrier_filt_state,
			mixed_rds_state,
			rds_demod_i_state,
			manchester_state,
			diff_state,
			frame_state,

			//MISC
			start_i_state,
			rds_pll
		);




		if (flag == true)
		{
			std::cerr << "Flag is One" << std::endl;
			if (queue.empty()) {
				std::cerr << "Queue is Empty" << std::endl;
				break;
			}
		}

		
	}

}



int main(int argc, char **argv)
{

	std::queue<std::shared_ptr<AudioBlock>> fm_queue;
	std::queue<std::shared_ptr<AudioBlock>> fm_queue_rds;
	std::mutex my_mutex;
	std::mutex my_mutex_rds;
	std::condition_variable my_cvar;
	std::condition_variable my_cvar_rds;
	Parameters &parameters = Parameters::getInstance();

	std::atomic<bool> stillReading(false);
	if (argc >= 2)
	{
		parameters.mode = atoi(argv[1]);
	}
	if (argc >= 3)
	{
		parameters.stereo = atoi(argv[2]);
	}
	if (argc >= 4)
	{
		parameters.logData = atoi(argv[3]);
	}

	std::cerr << "LOGGING = " << parameters.logData << std::endl;

	switch (parameters.mode)
	{
	case 0:
		parameters.rf_fs = 2400000.0;
		parameters.if_fs = 240000.0;
		parameters.audio_fs = 48000.0;

		parameters.rds = true;
		parameters.sps = 22; // sps * 2375?
		parameters.rdsUp = 209;
		parameters.rdsDown = 960;

		parameters.rf_down = 10;
		parameters.rf_up = 1;
		parameters.if_up = 1;
		parameters.if_down = 5;

		parameters.blockSize = BLOCK_SIZE * 100;		
		break;
	case 1:
		parameters.rf_fs = 1920000.0;
		parameters.if_fs = 240000.0;
		parameters.audio_fs = 48000.0;

		parameters.rf_down = 8;
		parameters.rf_up = 1;
		parameters.if_up = 1;
		parameters.if_down = 5;
		parameters.rds = false;

		parameters.blockSize = BLOCK_SIZE * 100;		
		break;
	case 2:
		parameters.rf_fs = 2400000.0;
		parameters.if_fs = 240000.0;
		parameters.audio_fs = 44100.0;

		parameters.rds = true;
		parameters.sps = 41; // sps * 2375?
		parameters.rdsUp = 209;
		parameters.rdsDown = 1920;

		parameters.rf_down = 10;
		parameters.rf_up = 1;
		parameters.if_up = 147;
		parameters.if_down = 800;

		parameters.blockSize = BLOCK_SIZE * 5 * 100;		
		break;
	case 3:
		parameters.rf_fs = 1152000.0;
		parameters.if_fs = 288000.0;
		parameters.audio_fs = 44100.0;

		parameters.rf_down = 4;
		parameters.rf_up = 1;
		parameters.if_up = 49;
		parameters.if_down = 320;
		parameters.rds = false;

		parameters.blockSize = BLOCK_SIZE * 1.5 * parameters.if_up;		

		break;
	default:
		std::cerr << "Invalid mode selected." << std::endl;
		std::cerr << "Mode " << parameters.mode << " has not been implemented" << std::endl;
		return -1;
	}

	parameters.carrier_lower_freq=18500;
	parameters.carrier_upper_freq=19500;
	parameters.stereo_lower_freq=22000;
	parameters.stereo_upper_freq=54000;

	std::cerr << "SDR Settings\n";
	std::cerr << "Mode: " << parameters.mode << std::endl;
	std::cerr << "rf_fs: " << parameters.rf_fs << std::endl;
	std::cerr << "if_fs: " << parameters.if_fs << std::endl;
	std::cerr << "audio_fs: " << parameters.audio_fs << std::endl;
	std::cerr << "rds: " << parameters.rds << std::endl;
	std::cerr << "sps: " << parameters.sps << std::endl;
	std::cerr << "stereo: " << parameters.stereo << std::endl;


	if(parameters.rds){
		std::thread rfThread = std::thread(rf_frontend, std::ref(parameters), std::ref(fm_queue),std::ref(fm_queue_rds), std::ref(my_mutex),std::ref(my_mutex_rds), std::ref(my_cvar),std::ref(my_cvar_rds), std::ref(stillReading));
		std::thread audioPaththread = std::thread(audio_path, std::ref(parameters), std::ref(fm_queue), std::ref(my_mutex), std::ref(my_cvar), std::ref(stillReading));
		std::thread rdbsPaththread = std::thread(rdbs_path, std::ref(parameters), std::ref(fm_queue_rds), std::ref(my_mutex_rds), std::ref(my_cvar_rds), std::ref(stillReading));
		rfThread.join();
		audioPaththread.join();
		rdbsPaththread.join();

	}else{
		std::thread rfThread = std::thread(rf_frontend, std::ref(parameters), std::ref(fm_queue),std::ref(fm_queue_rds), std::ref(my_mutex),std::ref(my_mutex_rds), std::ref(my_cvar),std::ref(my_cvar_rds), std::ref(stillReading));
		std::thread audioPaththread = std::thread(audio_path, std::ref(parameters), std::ref(fm_queue), std::ref(my_mutex), std::ref(my_cvar), std::ref(stillReading));
		rfThread.join();
		audioPaththread.join();
	}
	return 0;
}