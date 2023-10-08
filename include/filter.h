/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// add headers as needed
#include <iostream>
#include <vector>

// declaration of a function prototypes
void impulseResponseLPF(float, float, unsigned short int, std::vector<float> &);
void convolveFIR(std::vector<float> &, const std::vector<float> &, const std::vector<float> &);
void convolveAndResample(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, const int up, const int down);
void blockConvolution(const std::vector<float> &x, const std::vector<float> &IR, const int &m, std::vector<float> &state, std::vector<float> &y);
void convolveAndResampleBlock(const std::vector<float> &x, const std::vector<float> &IR, const int &m, std::vector<float> &state, std::vector<float> &y, const int up, const int down);
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
    const int down);
void fmDemod(std::vector<float> &fm_demod, const std::vector<float> &inPhase, float &inPhaseM1, const std::vector<float> &quadPhase, float &quadPhaseM1);

std::vector<float> my_lfilter(std::vector<float> &yb, const std::vector<float> &xb, const std::vector<float> &h, int block_size, std::vector<float> state);

void BPF(float fb, float fe,float fs, unsigned short int num_taps, std::vector<float> &h);

std::vector<float> impulseResponseRootRaisedCosine(int Fs, int N_taps);



#endif // DY4_FILTER_H
