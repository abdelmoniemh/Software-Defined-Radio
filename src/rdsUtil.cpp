#include <rdsUtil.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>

std::vector<std::vector<float>> P={
    {1,0,0,0,0,0,0,0,0,0},           //0
    {0,1,0,0,0,0,0,0,0,0},           //1
    {0,0,1,0,0,0,0,0,0,0},           //2
    {0,0,0,1,0,0,0,0,0,0},           //3
    {0,0,0,0,1,0,0,0,0,0},           //4
    {0,0,0,0,0,1,0,0,0,0},           //5
    {0,0,0,0,0,0,1,0,0,0},           //6
    {0,0,0,0,0,0,0,1,0,0},           //7
    {0,0,0,0,0,0,0,0,1,0},           //8
    {0,0,0,0,0,0,0,0,0,1},           //9
    {1,0,1,1,0,1,1,1,0,0},           //10
    {0,1,0,1,1,0,1,1,1,0},           //11
    {0,0,1,0,1,1,0,1,1,1},           //12
    {1,0,1,0,0,0,0,1,1,1},           //13
    {1,1,1,0,0,1,1,1,1,1},           //14
    {1,1,0,0,0,1,0,0,1,1},           //15
    {1,1,0,1,0,1,0,1,0,1},           //16
    {1,1,0,1,1,1,0,1,1,0},           //17
    {0,1,1,0,1,1,1,0,1,1},           //18
    {1,0,0,0,0,0,0,0,0,1},           //19
    {1,1,1,1,0,1,1,1,0,0},           //20
    {0,1,1,1,1,0,1,1,1,0},           //21
    {0,0,1,1,1,1,0,1,1,1},           //22
    {1,0,1,0,1,0,0,1,1,1},           //23
    {1,1,1,0,0,0,1,1,1,1},           //24
    {1,1,0,0,0,1,1,0,1,1}            //25
};

std::unordered_map<std::string, std::vector<int>> Syndromes = {
    {"A", {1,1,1,1,0,1,1,0,0,0}},
    {"B", {1,1,1,1,0,1,0,1,0,0}},
    {"C", {1,0,0,1,0,1,1,1,0,0}},
    {"C\'", {1,1,1,1,0,0,1,1,0,0}},
    {"D", {1,0,0,1,0,1,1,0,0,0}},
};

void manchesterDecoding(
    const std::vector<float> &toBeDecoded,
    std::vector<bool> &cdr,
    std::vector<bool> &state)
{

    int i = 0;
    cdr = std::vector<bool>();

    while (i < toBeDecoded.size() - 1)
    {

        if ((toBeDecoded[i] > 0 && toBeDecoded[i + 1] > 0) || (toBeDecoded[i] < 0 && toBeDecoded[i + 1] < 0))
        {
            cdr.push_back(false);
            i += 2;
            continue;
        }
        if (i + 1 >= toBeDecoded.size())
        {
            state = std::vector<bool>(1, toBeDecoded.back());
            break;
        }
        if (toBeDecoded[i] < toBeDecoded[i + 1])
        {
            cdr.push_back(false);
            i += 2;
        }
        else
        {
            cdr.push_back(true);
            i += 2;
        }
    }

    state = std::vector<bool>();
}

void differentialDecoding(
    const std::vector<bool> &toBeDecoded,
    std::vector<bool> &postCdr,
    std::vector<bool> &state)
{
    postCdr = std::vector<bool>();
    std::vector<bool> toDecode = state;
    toDecode.insert(
        toDecode.end(),
        toBeDecoded.begin(),
        toBeDecoded.end()
    );

    for (int i = 0; i<toDecode.size()-2; i++){
        postCdr.push_back(toDecode[i]^toDecode[i+1]);
    }
    state = std::vector<bool>(1, toDecode.back());
}

void frameSync(
    const std::vector<int> &cdrOut,
    std::vector<int> &state)
{
    std::vector<int> block = std::vector<int>();
    block.reserve(104);
    for (int x : state){
        block.push_back(x);
    }
    for (int i = 0; i<cdrOut.size(); i++){
        block.push_back(cdrOut[i]);
        if (block.size() == 104){

            auto start=block.begin();

            //check all 26 bit blocks
            std::vector<int> syndrome_A = getSyndrome(std::vector<int>(start, start + 26));       //Check A
            start+=26;
            std::vector<int> syndrome_B = getSyndrome(std::vector<int>(start, start + 26));       //Check B
            start+=26;
            std::vector<int> syndrome_C = getSyndrome(std::vector<int>(start, start + 26));       //Check C/C'
            start+=26;
            std::vector<int> syndrome_D = getSyndrome(std::vector<int>(start, start + 26));       //Check D


            std::string offsetType_A = matchSyndrome(syndrome_A);
            std::string offsetType_B = matchSyndrome(syndrome_B);
            std::string offsetType_C = matchSyndrome(syndrome_C);
            std::string offsetType_D = matchSyndrome(syndrome_D);



            // if (offsetType_A.compare("?") || offsetType_B.compare("?") || offsetType_C.compare("?") || offsetType_D.compare("?")){
            //     std::cerr << "Matched block " << offsetType_A<<" "<<offsetType_B <<" "<<offsetType_C<<" "<<offsetType_D<< std::endl;
            // }

            int PI=0;
            

            // DATA DECODING
            //compare returns 0 if valid
            if(!offsetType_A.compare("A")){  // match for A
                PI=arrayToInt(std::vector<int>(block.begin(), block.begin() + 16));
                std::cerr<<std::hex<<"Block A MATCH. PI CODE: "<<PI<< std::dec <<std::endl;
            }

            int PTY=0;
            int GTY=0;

            if(!offsetType_B.compare("B")){  // match for B
                PTY=arrayToInt(std::vector<int>(block.begin()+26, block.begin() + 26 + 4));
                GTY=arrayToInt(std::vector<int>(block.begin()+26+7, block.begin() + 26+7+5));

                std::cerr<<std::hex<<"Block B MATCH. PTY CODE: "<<PTY<<std::dec <<std::endl;
                std::cerr<<std::hex<<"Block B MATCH. GTY CODE: "<<GTY<<std::dec <<std::endl;
            }

            if(!offsetType_C.compare("C") || !offsetType_C.compare("C\'")){  // match for C
                // Irrelevant control bits...
                std::cerr << "Block C/C' MATCH." << std::endl;
            }

            if(!offsetType_D.compare("D")){  // match for D
                std::cerr << "Block D MATCH." << std::endl;
            }


            if(offsetType_A=="A" &&
               offsetType_B=="B" &&
               (offsetType_C=="C" || offsetType_A=="C\'")   &&
               offsetType_A=="D"
            ){
                std::cerr << "Matched Frame " << offsetType_A<<" "<<offsetType_B <<" "<<offsetType_C<<" "<<offsetType_D<< std::endl;
                block.clear();
            }
            else{
                block.erase(block.begin());
            }


        }
    }
    state = block;
}

std::vector<int> getSyndrome(std::vector<int> block) {
    std::vector<int> matrix_col(26);
    std::vector<int> result(10);

    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 26; j++) {
            matrix_col[j] = P[j][i];
        }
        
        std::vector<bool> block_bool(block.begin(), block.end());
        std::vector<bool> matrix_col_bool(matrix_col.begin(), matrix_col.end());
        
        std::vector<bool> prod(26);
        transform(block_bool.begin(), block_bool.end(), matrix_col_bool.begin(), prod.begin(), 
            [](bool a, bool b){ return a & b; });

        for (bool elem : prod) {
            result[i] ^= elem;
        }
    }
    return result;
}

std::string matchSyndrome(std::vector<int> syndrome){
    std::vector<std::string> blocks{"A", "B", "C", "C\'", "D"};
    for(auto elem:blocks){
        if(syndrome==Syndromes[elem]){
            return elem;
        }
    }
    return "?";
}



int findLocalMaxMin(std::vector<float> data, int symbol_rate){
    int max = -100000; int min = 100000;
    int max_i = 0; int min_i = 0;
    for (int i = 0; i < symbol_rate; i++){
        if (data[i] > max){
            max_i = i; max = data[i];
        }
        if (data[i] < min){
            min_i = i; min = data[i];
        }
    }
    if (abs(max) > abs(min)){
        return max_i;
    } else {
        return min_i;
    }
}



int arrayToInt(std::vector<int> arr){
    int result=0;
    for(int i=arr.size()-1;i>=0;i--){
        result+=arr[i]*pow(2,i);
    }
    return result;

}