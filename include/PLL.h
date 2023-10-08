#ifndef PLL_H_
#define PLL_H_

#include <vector>
using namespace std;

// NOTE: since PLL is handeled in a a class, state saving is handeled "automatically"

class PLL{

    std::vector<float> ncoOut;

    // The members of PLL are the saved state    

    float s_integrator;
    float s_phaseEst;
    float s_feedback_I;
    float s_feedback_Q;
    float s_ncoOut;
    float s_trigOffset;


public:

    PLL(int blockSize);

    //TODO: change pllIn to block object
    std::vector<float> lock_phase_block(const vector<float> &pllIn, const float &freq, const float &Fs, 
    const float &ncoScale, const float &phaseAdjust, const float &normBandwidth);    

};

#endif