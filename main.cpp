
#include "indexes.hpp"
#include <sdsl/testutils.hpp>
#include "ED.hpp"
#include <ConfCheckFMI.hpp>
#include <algorithm>
#include <chrono>
int main(){
    /*int requiredArguments = 2; // index and patterns

    if (argc < requiredArguments) {
        std::cerr << "Insufficient number of arguments" << std::endl;
        //showUsage();
        return EXIT_FAILURE;
    }
    int maxED = 0;

    for (int i = 1; i < argc - requiredArguments; i++){

    }*/
    std::string log_file_name = "Resources/Unique/uniqueAO.txt";
    std::string model_file_name = "Resources/Model/AOModel.txt";
    //std::string log_file_name = "Resources/logs_str.txt";
    //std::string model_file_name = "Resources/model_str.txt";
    //std::string log_file_name = "Resources/traces_cf_B.txt";
    //std::string model_file_name = "Resources/B_Model_cases.txt";
    ConfCheckFMI ccf;

    auto index_str = ccf.load_model(model_file_name);
    auto patterns = ccf.load_log(log_file_name);
    const unsigned char * FMIndex_str= reinterpret_cast <const unsigned char *> (index_str.c_str());

    auto main_index = index_bidirectional_waveletindex(FMIndex_str);
    int k = 1;
    auto as = ED(&main_index, FMIndex_str);
    bool fw = true, complete = true, all=true, replace = true;

    std::vector<int> foundPerK(k + 2);
    std::chrono::high_resolution_clock::time_point start, end;
    int totalTime = 0;
    
    for(const auto &pattern : patterns){
        start = std::chrono::high_resolution_clock::now();  
        int pK = as.CheckConf(pattern, k,fw, complete, replace);
        end = std::chrono::high_resolution_clock::now();
        totalTime += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        if (pK == 1){
            std::cout << pattern <<std::endl;
            // as.printOcc(fw == hamming, pattern);
        //    std::cout << "Not found:" << pattern << std::endl;
        }
        foundPerK[pK != -1 ? pK : k + 1]++;

        

    }
    for (int i = 0; i <= k; i++){
        std::cout<< i <<": "<< foundPerK[i] << std::endl;
    }
    std::cout<< "Not found: "<< foundPerK[k + 1] << std::endl;
    std::cout<< "Total number of patterns: "<< patterns.size() << std::endl;
    std::cout<< "Total time (ms): "<< totalTime << std::endl;
    std::cout<< "Average time per pattern (ms): "<< totalTime/patterns.size()<< std::endl;
    std::cout<< "Total time (s): "<< totalTime/(1000.0)<< std::endl;
}