#include <vector>
#include <iostream>
#include <cmath>

int arrayToInt(std::vector<int> arr){
    int result=0;
    for(int i=arr.size()-1;i>=0;i--){
        result+=arr[i]*pow(2,i);
    }
    return result;

}


int main(){

    std::vector<int> a{1,1,1};

    std::cerr<<std::hex<<arrayToInt(a)<<"\n";


}