#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "audioBlock.h"
#include "PLL.h"

#define BLOCK_SIZE 1024

void processBlock(std::vector<float> &block_data);

