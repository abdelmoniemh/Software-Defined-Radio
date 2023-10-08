#include <global.h>
#include <audioBlock.h>
#include "estimatePSD.cpp"
#include <PLL.h>
#include <complex>

#include <rdsUtil.h>
#include <algorithm>

Parameters &parameters = Parameters::getInstance();

AudioBlock::AudioBlock(std::vector<float> &data, int id, int size)
{
    blockData = data;
    blockId = id;
    blockSize = size;
}

void AudioBlock::preprocessRFData(
    std::vector<float> &stateI,
    std::vector<float> &stateQ,
    float &inPhaseM1,
    float &quadPhaseM1,
    const std::vector<float> &rf_coeff,
    const int rf_up,
    const int rf_down)
{
    std::vector<float> i_data, q_data;
    split_audio_into_channels(blockData, i_data, q_data);
    std::vector<float> i_filt(blockSize / 2, 0.0), q_filt(blockSize / 2, 0.0);

    multiConvolveAndResampleBlock(
        i_data,
        q_data,
        rf_coeff,
        rf_coeff,
        blockId,
        stateI,
        stateQ,
        i_filt,
        q_filt,
        rf_up,
        rf_down);

    if (parameters.logData)
    {
        estimatePSD(i_filt, parameters.if_fs, "../data/PSD/mode_" + std::to_string(parameters.mode) + "/i_data/" + "block_" + std::to_string(blockId) + +".txt");
        estimatePSD(q_filt, parameters.if_fs, "../data/PSD/mode_" + std::to_string(parameters.mode) + "/q_data/" + "block_" + std::to_string(blockId) + +".txt");
    }

    fmDemod(fm_demod, i_filt, inPhaseM1, q_filt, quadPhaseM1);

    if (parameters.logData)
    {
        estimatePSD(fm_demod, parameters.if_fs, "../data/PSD/mode_" + std::to_string(parameters.mode) + "/fm_demod/" + "block_" + std::to_string(blockId) + +".txt");
    }
    std::cerr << "Processed rf data for block " << blockId << std::endl;
}

void AudioBlock::processMonoAudio(
    std::vector<float> &audio_coeff,
    std::vector<float> &delay_coeffs,
    std::vector<float> &filterState,
    std::vector<float> &stateDelay,
    std::queue<float> &monoQ,
    const int if_up,
    const int if_down)
{

    monoData = std::vector<float>(fm_demod.size(), 0.0);
    std::vector<float> fmDemodDelayed(fm_demod.size(), 0.0);

    if (parameters.stereo)
    {

        for (float fmElem : fm_demod)
        {
            monoQ.push(fmElem);
        }
        monoData.clear();
        monoData.resize(fm_demod.size(), 0.0);
        for (int i = 0; i < fm_demod.size(); i++)
        {
            if (!monoQ.empty())
            {
                fmDemodDelayed[i] = monoQ.front();
                monoQ.pop();
            }
        }

        convolveAndResampleBlock(fmDemodDelayed, audio_coeff, blockId, filterState, monoData, if_up, if_down);
    }
    else
    {
        std::cerr << "Size of fm_demod (num of taps cant exceed this) " << fm_demod.size() << std::endl;
        convolveAndResampleBlock(fm_demod, audio_coeff, blockId, filterState, monoData, if_up, if_down);
    }

    // if (parameters.logData)
    // {
    //     estimatePSD(monoData, parameters.audio_fs, "../data/PSD/mode_" + std::to_string(parameters.mode) + "/monodata/" + "block_" + std::to_string(blockId) + ".txt");
    // }

    if (blockId == 0)
    {
        logforPSD("_monoData0", monoData);
        logforPSD("_fm_demod0", fm_demod);
    }

    if (blockId == 1)
    {
        logforPSD("_monoData1", monoData);
        logforPSD("_fm_demod1", fm_demod);
    }

    if (blockId == 50)
    {
        logforPSD("_monoData50", monoData);
        logforPSD("_fm_demod50", fm_demod);
    }

    std::cerr << "Processed mono data for block " << blockId << std::endl;
}

void AudioBlock::processStereoAudio(
    std::vector<float> &audio_coeffs,
    std::vector<float> &pilot_coeff,
    std::vector<float> &stereo_coeff,

    // States
    std::vector<float> &stateCarrier,
    std::vector<float> &stateStereo,
    std::vector<float> &stateStereoMixed,

    const int if_up,
    const int if_down,
    const int if_fs,

    // PLL
    PLL &pll)
{

    // Exract Pilot
    std::vector<float> pilotData(fm_demod.size(), 0.0);
    std::vector<float> stereoChannel(fm_demod.size(), 0.0);
    std::vector<float> pllOut(fm_demod.size(), 0.0);
    std::vector<float> mixed_signal(fm_demod.size(), 0.0);
    std::vector<float> mixed_signal_filtered(fm_demod.size() / if_down * if_up, 0.0);

    multiConvolveAndResampleBlock(
        fm_demod,     // input 1
        fm_demod,     // input 2
        pilot_coeff,  // coeff 1
        stereo_coeff, // coeff 2
        blockId,
        stateCarrier,  // state 1
        stateStereo,   // state 2
        pilotData,     // output 1
        stereoChannel, // output 2
        1,
        1);

    // Lock pilot phase (PLL)
    pllOut = pll.lock_phase_block(pilotData, 19000, 240000, 2.0, 0, 0.01);

    // Mix
    for (int i = 0; i < stereoChannel.size(); i++)
    {
        mixed_signal[i] = pllOut[i] * stereoChannel[i] * 2;
    }

    convolveAndResampleBlock(mixed_signal, audio_coeffs, blockId, stateStereoMixed, stereoData, if_up, if_down);

    if (this->blockId == 0)
    {
        logforPSD("_fm_demod0", fm_demod);
        logforPSD("_pilot_coeffs0", pilot_coeff);
        logforPSD("_stereo_coeffs0", stereo_coeff);
        logforPSD("_pilotData0", pilotData);
        logforPSD("_stereoChannel0", stereoData);
        logforPSD("_pllOut0", pllOut);
        logforPSD("_mixed_signal0", mixed_signal);

        logforPSD("_stateCarrier0", stateCarrier);
    }

    if (this->blockId == 1)
    {
        logforPSD("_fm_demod1", fm_demod);
        logforPSD("_pilot_coeffs1", pilot_coeff);
        logforPSD("_stereo_coeffs1", stereo_coeff);
        logforPSD("_pilotData1", pilotData);
        logforPSD("_stereoChannel1", stereoData);
        logforPSD("_pllOut1", pllOut);
        logforPSD("_mixed_signal1", mixed_signal);
    }

    std::cerr << "Processed stereo data for block " << blockId << std::endl;
}

void AudioBlock::processRBDS(
    //Filters
    std::vector<float> &rds_coeff,
    std::vector<float> &rds_carrier_coeffs,
    std::vector<float> &rds_demod_coeffs,
    std::vector<float> &cosine_coeffs,

    //DELAY QUEUE
    std::queue<float> &rds_Q,

    //FILTER STATES
    std::vector<float> &rds_filt_state,
    std::vector<float> &rds_carrier_filt_state,
    std::vector<float> &mixed_rds_state,
    std::vector<float> &rds_demod_i_state,
    std::vector<bool> &manchester_state,
    std::vector<bool> &diff_state,
    std::vector<int> &frame_state,

    //MISC
    int &start_i_state,
    PLL &pll)
{


    //std::cerr << "Processing RBDS for block "<< blockId<<"\n";



    std::vector<float> rds_filt;
    blockConvolution(fm_demod, rds_coeff, blockId, rds_filt_state, rds_filt); // extract rds channel

    for (float x : rds_filt)
    {
        rds_Q.push(x);
    }

    rds_filt = std::vector<float>(fm_demod.size(), 0.0);
    for (int i = 0; i < fm_demod.size(); i++)
    {
        rds_filt[i] = rds_Q.front();
        rds_Q.pop();
    }

    std::vector<float> rds_carrier(fm_demod.size(), 0.0);
    for (int i = 0; i < rds_filt.size(); i++)
    {
        rds_carrier[i] = rds_filt[i] * rds_filt[i]; // Square signal
    }

    std::vector<float> rds_carrier_filt(fm_demod.size(),0.0);
    blockConvolution(rds_carrier, rds_carrier_coeffs, blockId, rds_carrier_filt_state, rds_carrier_filt); // extract rds carrier wave

    std::vector<float> pllOut;
    pllOut = pll.lock_phase_block(rds_carrier_filt, 114e3, 240000, 0.5, 0, 0.01); // this is rds_carrier_i



    std::vector<float> rds_mixed_i(fm_demod.size(), 0.0);
    for (int i = 0; i < rds_filt.size(); i++)
    {
        rds_mixed_i[i] = 2.0 * pllOut[i] * rds_filt[i];
    }

    std::vector<float> rds_mixed_filt_i;
    convolveAndResampleBlock(rds_mixed_i, rds_demod_coeffs, blockId, mixed_rds_state, rds_mixed_filt_i, parameters.rdsUp, parameters.rdsDown); // extract rds channel


    std::vector<float> rds_demod_i;
    blockConvolution(rds_mixed_filt_i, cosine_coeffs, blockId, rds_demod_i_state, rds_demod_i); // extract rds channel


    // DATA EXTRACTION AND CLOCK RECOVERY

    int delay = (int)(((151.0) * (float)parameters.rdsUp / (float)parameters.rdsDown) + (151.0 - 1.0) / 2.0);
    rds_demod_i = std::vector<float>(rds_demod_i.begin() + delay, rds_demod_i.end());   // REMOVE DELAY DUE TO RRC FILTER

    int start_i;
    if (blockId == 0)
    {
        start_i = findLocalMaxMin(rds_demod_i, parameters.sps); // FIND FIRST SAMPLE POINT FOR SPS SKIPPNIG 
    }
    else
    {
        start_i = start_i_state;
    }

    rds_demod_i = std::vector<float>(rds_demod_i.begin() + start_i, rds_demod_i.end());

    std::vector<float> toBeDecoded = std::vector<float>();
    toBeDecoded.reserve(rds_demod_i.size()/(float)parameters.sps);
    int index = -1;

    while (index + 1 < rds_demod_i.size())
    {
        index++;
        if (index % 11 == 0)
        {
            toBeDecoded.push_back(rds_demod_i[index]);
        }
    }
    start_i_state = 11 - index % 11;

    std::vector<bool> cdr;
    manchesterDecoding(toBeDecoded, cdr, manchester_state);

    std::vector<bool> postCdr;
    differentialDecoding(cdr, postCdr, diff_state);

    std::vector<int> frameInput= std::vector<int>();
    frameInput.reserve(postCdr.size());
    for (bool x : postCdr){
        frameInput.push_back((int)x);
    }

    frame_state.reserve(104);
    frameSync(frameInput, frame_state);



    

//    std::cerr << "Processed RBDS for block " << blockId << std::endl;

    // for block in stream
        // do matrix multiply
            // if match, print block and data and empty block



}

void AudioBlock::outputAudio(
    bool stereo,
    bool toTerminal,
    std::vector<float> *allDataLeft,
    std::vector<float> *allDataRight)
{
    std::vector<float> toOutput(monoData.size(), 0.0);
    std::vector<float> left_channel(monoData.size(), 0.0);
    std::vector<float> right_channel(monoData.size(), 0.0);

    if (stereo)
    {
        for (unsigned int k = 0; k < monoData.size(); k++)
        {
            // stereoData[k]*=32767/2;
            toOutput[k] = monoData[k] * 32767 / 2;
            left_channel[k] = (monoData[k] + stereoData[k]) * 32767.0 / 2.0;
            right_channel[k] = (monoData[k] - stereoData[k]) * 32767.0 / 2.0;
        }
        // toOutput=left_channel;
    }
    else
    { // Only mono data
        for (unsigned int k = 0; k < monoData.size(); k++)
        {
            toOutput[k] = monoData[k] * 32767 / 2;
        }
    }

    if (toTerminal)
    {
        if (!stereo)
        {
            short int sample;
            for (unsigned int k = 0; k < toOutput.size(); k++)
            {
                if (std::isnan(toOutput[k]))
                    sample = 0;
                else
                    sample = static_cast<short int>(toOutput[k]);
                fwrite(&sample, sizeof(short int), 1, stdout);
            }
        }
        else
        {
            short int sample;
            for (unsigned int k = 0; k < toOutput.size() * 2; k++)
            {
                float toOutput = k % 2 == 0 ? left_channel[k / 2] : right_channel[k / 2];
                if (std::isnan(toOutput))
                    sample = 0;
                else
                    sample = static_cast<short int>(toOutput);
                fwrite(&sample, sizeof(short int), 1, stdout);
            }
        }
    }

    if (stereo && allDataLeft && allDataRight)
    {
        for (unsigned int k = 0; k < toOutput.size(); k++)
        {
            // allDataLeft->push_back(stereoData[k] * 32767 / 2);
            // allDataRight->push_back(stereoData[k] * 32767 / 2);

            allDataLeft->push_back(left_channel[k]);
            allDataRight->push_back(right_channel[k]);
        }
    }
    else if (allDataLeft)
    {
        for (unsigned int k = 0; k < toOutput.size(); k++)
        {
            allDataLeft->push_back(toOutput[k]);
        }
    }
}

int AudioBlock::getBlockId()
{
    return blockId;
}