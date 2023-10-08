#ifndef RDSUTILS_H_
#define RDSUTILS_H_

#include <vector>
#include <unordered_map>
#include <bitset>


void manchesterDecoding(
    const std::vector<float>& toBeDecoded,
    std::vector<bool>& cdr,
    std::vector<bool>& state);

void differentialDecoding(
    const std::vector<bool>& toBeDecoded,
    std::vector<bool>& postCdr,
    std::vector<bool>& state);


std::vector<int> getSyndrome(std::vector<int> block);
std::string matchSyndrome(std::vector<int> syndrome);



void frameSync(
    const std::vector<int>& cdrOut,
    std::vector<int>& state);

int findLocalMaxMin(std::vector<float> data, int symbol_rate);

int arrayToInt(std::vector<int> arr);

#endif