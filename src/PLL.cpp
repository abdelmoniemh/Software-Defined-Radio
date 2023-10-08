#include "PLL.h"
#include "dy4.h"
#include <cmath>

#include <iostream>


PLL::PLL(int blockSize){

    //std::vector<float> pllIn(blockSize, 0.0);   // change this to a block object

    ncoOut = std::vector<float>(blockSize+1, 0.0);
    ncoOut[0]=1.0;

    // The members of PLL are the saved state

    float s_integrator=0.0;
    float s_phaseEst=0.0;
    float s_feedback_I=1.0;
    float s_feedback_Q=0.0;
    float s_ncoOut=1.0;
    float s_trigOffset=0.0;


}


std::vector<float> PLL::lock_phase_block(const vector<float> &pllIn, const float &freq, const float &Fs, 
                           const float &ncoScale, const float &phaseAdjust, const float &normBandwidth){

    float Cp = 2.666;
    float Ci = 3.555;

    float Kp = normBandwidth*Cp;

    float Ki = (normBandwidth*normBandwidth)*Ci;


    float error_I=0.0;
    float error_Q=0.0;

    int max=pllIn.size();

    ncoOut[0]=s_ncoOut;

    for(int k=0;k<pllIn.size();k++){
        error_I=pllIn[k]*(+s_feedback_I);
        error_Q=pllIn[k]*(+s_feedback_Q);

        float error_D=std::atan2(error_I,error_Q);

        s_integrator+=Ki*error_D;

        s_phaseEst+=+Kp*error_D+s_integrator;

        s_trigOffset+=1;
        float trigArg=2*PI*(freq/Fs)*(s_trigOffset)+s_phaseEst;

        s_feedback_I=std::cos(trigArg);
        s_feedback_Q=std::sin(trigArg);

        ncoOut[k+1]=std::cos(trigArg*ncoScale+phaseAdjust);
    }


    s_ncoOut=ncoOut[pllIn.size()-1];

    return ncoOut;

}

