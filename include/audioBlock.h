#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "PLL.h"
#include <queue>

class AudioBlock {
    private:
        std::vector<float> blockData;
        int blockId;
        int blockSize;

        std::vector<float> fm_demod;
        std::vector<float> monoData;
        std::vector<float> stereoData;
        std::vector<float> left_channel;
        std::vector<float> right_channel;

        

    public:
        AudioBlock(std::vector<float>& data, int blockId, int size);


        void preprocessRFData(
            std::vector<float>& stateI,
            std::vector<float>& stateQ,
            float& inPhaseM1,
	        float& quadPhaseM1,
            const std::vector<float>& rf_coeff,
            const int rf_up,
            const int rf_down
        );

        void processMonoAudio(
            std::vector<float>& audio_coeff,
            std::vector<float>& delay_coeffs,
            std::vector<float>& filterState,
            std::vector<float>& stateDelay,
            std::queue<float>& monoQ,

            const int if_up,
            const int if_down
            );

        void processStereoAudio(
            std::vector<float>& audio_coeffs,
            std::vector<float>& pilot_coeff,
            std::vector<float>& stereo_coeff,

            //States
            std::vector<float>& stateCarrier,
            std::vector<float>& stateStereo,
            std::vector<float>& stateStereoMixed,

            const int if_up,
            const int if_down,
            const int if_fs,

            //PLL
            PLL &pll
        );

        void processRBDS(    
            std::vector<float> &rds_coeff,
            std::vector<float> &rds_carrier_coeffs,
            std::vector<float> &rds_demod_coeffs,
            std::vector<float> &cosine_coeffs,
            std::queue<float>& rds_Q,
            std::vector<float> &rds_filt_state,
            std::vector<float> &rds_carrier_filt_state,
            std::vector<float> &mixed_rds_state,
            std::vector<float> &rds_demod_i_state,
            std::vector<bool> &manchester_state,
            std::vector<bool> &diff_state,
            std::vector<int> &frame_state,
            int& start_i_state,
            PLL &pll);
            


        void outputAudio(
            bool stereo,
            bool toTerminal,
            std::vector<float>* allDataLeft = nullptr,
            std::vector<float>* allDataRight = nullptr
        );

        int getBlockId();

};
