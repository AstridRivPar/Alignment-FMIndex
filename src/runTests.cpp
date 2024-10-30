
#include "indexes.hpp"
#include <sdsl/testutils.hpp>
#include "ED.hpp"
#include <ConfCheckFMI.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <boost/algorithm/string/join.hpp>

Index makeIndex(const unsigned char * FMIndex_str){    
    auto main_index = index_bidirectional_waveletindex(FMIndex_str);
    return main_index;
}

void Log2LogOptimalStats(const std::vector<std::string> &traces, int maxK, std::ofstream &of){
   
    
    std::vector<int> foundPerK(maxK + 2);
    // of << "T No.\tTrace Length \t Optimal Value \t No. of Opt" << std::endl;
    

    for(size_t i = 0; i < traces.size(); i++){
        std::cout << i <<std::endl;
        //build index
        std::stringstream ss;
        ss << ",";
        for (size_t j = 0; j < traces.size(); j++){
            if(j != i){
                ss << traces[j] <<",";
            }
        }
        std::string index_str = ss.str();
        
        const unsigned char * FMIndex_str= reinterpret_cast <const unsigned char *> (index_str.c_str());
        auto main_index = index_bidirectional_waveletindex(FMIndex_str);
        ED ed = ED(&main_index, FMIndex_str);
        int pK = ed.CheckConf(traces[i],maxK,true,true,false);
        
        of << i + 1<<"\t";
        of << traces[i].size() << "\t";
        if (pK != -1){
            of << pK << "\t";
            of << ed.getNoOpt(traces[i]) <<"\t";
        }
        ed.getBest(of);
        of << std::endl;
        foundPerK[pK != -1 ? pK : maxK + 1]++;
        
           
    } 
    
    for (int i = 0; i <= maxK; i++){
        of << i <<":"<< foundPerK[i] << std::endl;
    }
    of<< "NF:" <<foundPerK[maxK + 1] <<std::endl;     

}
double stdev(const std::vector<double>& vals) {
  double sum = std::reduce(vals.begin(), vals.end());
  double mean = sum / static_cast<double>(vals.size());

  // variance
  double squaredDifference = 0.0;
  for (unsigned int i = 0; i < vals.size(); i++)
    squaredDifference += std::pow(vals[i] - mean, 2);
  

  return std::sqrt(squaredDifference / static_cast<double>(vals.size()));
}

void OneTest(bool replace, bool fw, int k, ED &ed, const std::vector<std::string> &traces,  std::vector<int> &foundPerK){

    for(const auto &trace : traces){
        int pK = ed.CheckConf(trace, k, fw, true, replace);
        foundPerK[pK != -1 ? pK : k + 1]++;            
    }
}
void DoTests(bool replace, bool fw, bool complete, int maxED, Index &idx, const std::vector<std::string> &traces, const unsigned char * FMIndex_str, std::ostream &of, const int numberOfTests = 5){
    std::chrono::high_resolution_clock::time_point start, end;
    std::vector<double> meanV(maxED + 1);
    std::vector<double> stdevV(maxED + 1);
    //auto k = maxK;
    for (int k = 0; k <= maxED; k++){
        std::cout << "K = " << k << std::endl;
        std::vector<double> times(numberOfTests);
        std::vector<int> foundPerK(k + 2);
        
        of << "K:" << k << std::endl;
        for(int i = 0; i < numberOfTests; i++){
            std::cout << "Test " << i + 1 << std::endl;
            std::fill(foundPerK.begin(), foundPerK.end(), 0);
            ED ed = ED(&idx, FMIndex_str);
            start = std::chrono::high_resolution_clock::now();
            OneTest(replace, fw, k, ed, traces, foundPerK);
            end = std::chrono::high_resolution_clock::now();
            int64_t totalTime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            // std::cout << totalTime <<std::endl;
            times[i] = (totalTime) /(1000.0);  
            
           

        }
        sort(times.begin(), times.end());
        //erase highest and lowest time
        times.pop_back();
        times.erase(times.begin());


        //calculate mean
  
        double mean = std::reduce(times.begin(), times.end())/static_cast<double>(times.size());
        meanV[k] = mean;

        //calculate standard deviation
        double standardDev = stdev(times);
        stdevV[k] = standardDev;

        //print header
        for (int i = 0; i < foundPerK.size() - 1; i++){
            of << i << "\t";
        }
        of << "NF\t"<< std::endl;

        //print number of traces with k mismatches
        for (int i = 0; i <= k; i++){
            of << foundPerK[i] << "\t";
        }
        of<< foundPerK[k + 1] << "\t" <<std::endl;


        //print times
        of << "Times (s)" << std::endl;
        for (int i = 0; i < times.size(); i++){
            of << times[i]<<"\t" ;
        }
        of << std::endl;

        
        of << "Mean: " << mean << std::endl;
        of << "SD: " << standardDev << std::endl;
        of << "------------" <<std::endl;
    }



    //print number of traces with k mismatches
    of << "Summary of Times (s) \nMean \t STDev\n"; 
    for (int i = 0; i <= maxED; i++){
        of << meanV[i] << "\t" << stdevV[i] <<std::endl;
    }
    
    of << "************************************" << std::endl;
    of << std::endl;
}
void DoTestsL2L(bool replace, bool fw, bool complete, int maxED, const std::vector<std::string> &traces, std::ostream &of, const int numberOfTests = 5){
    std::chrono::high_resolution_clock::time_point start, end;
    std::vector<double> meanV(maxED + 1);
    std::vector<double> stdevV(maxED + 1);
    int64_t totalTime = 0;
    //auto k = maxK;
    for (int k = maxED; k <= maxED; k++){
        std::cout << "K = " << k << std::endl;
        std::vector<double> times(numberOfTests);
        std::vector<int> foundPerK(k + 2);
        
        of << "K:" << k << std::endl;
        for(int i = 0; i < numberOfTests; i++){
            totalTime = 0;
            std::cout << "Test " << i + 1 << std::endl;
            std::fill(foundPerK.begin(), foundPerK.end(), 0);
           
            
            //
            
            
            for(size_t j = 0; j < traces.size(); j++){
                //build index
                std::stringstream ss;
                ss << ",";
                for (size_t k = 0; k < traces.size(); k++){
                    if(k != j){
                        ss << traces[k] <<",";
                    }
                }
                std::string index_str = ss.str();
                // std::vector<std::string> idx(traces.begin(), traces.end());
                // // idx.erase(idx.begin()+j);
                // std::string index_str = ","+traces[j]+",";
                // if (j==36) 
                //     std::cout << index_str << std::endl;
                // std::string index_str = ","+boost::algorithm::join(traces, ",");
                const unsigned char *FMIndex_str= reinterpret_cast <const unsigned char *> (index_str.c_str());
                Index main_index = index_bidirectional_waveletindex(FMIndex_str);
                start = std::chrono::high_resolution_clock::now();
                ED ed(&main_index, FMIndex_str);
                int pK = ed.CheckConf(traces[j],maxED,fw,complete,replace);
                end = std::chrono::high_resolution_clock::now();
                totalTime += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
                foundPerK[pK != -1 ? pK : k + 1]++;

            }
            
            // std::cout << totalTime <<std::endl;
            times[i] = (totalTime) /(1000.0); 

        }
        sort(times.begin(), times.end());
        //erase highest and lowest time
        times.pop_back();
        times.erase(times.begin());


        //calculate mean
  
        double mean = std::reduce(times.begin(), times.end())/static_cast<double>(times.size());
        meanV[k] = mean;

        //calculate standard deviation
        double standardDev = stdev(times);
        stdevV[k] = standardDev;

        //print header
        for (int i = 0; i < foundPerK.size() - 1; i++){
            of << i << "\t";
        }
        of << "NF\t"<< std::endl;

        //print number of traces with k mismatches
        for (int i = 0; i <= k; i++){
            of << foundPerK[i] << "\t";
        }
        of<< foundPerK[k + 1] << "\t" <<std::endl;


        //print times
        of << "Times (s)" << std::endl;
        for (int i = 0; i < times.size(); i++){
            of << times[i]<<"\t" ;
        }
        of << std::endl;

        
        of << "Mean: " << mean << std::endl;
        of << "SD: " << standardDev << std::endl;
        of << "------------" <<std::endl;
    }



    //print number of traces with k mismatches
    of << "Summary of Times (s) \nMean \t STDev\n"; 
    for (int i = 0; i <= maxED; i++){
        of << meanV[i] << "\t" << stdevV[i] <<std::endl;
    }
    
    of << "************************************" << std::endl;
    of << std::endl;
}
/**
 * Test Indexes
 */
void TestIndexes(std::string index_str, std::ofstream &of){
    const unsigned char * FMIndex_str= reinterpret_cast <const unsigned char *> (index_str.c_str());
    
    std::chrono::high_resolution_clock::time_point start, end;
    std::vector<double> times(5);
    for (int i = 0;i < times.size(); i++){
        start = std::chrono::high_resolution_clock::now();
        Index idx = makeIndex(FMIndex_str);
        end = std::chrono::high_resolution_clock::now();
        int64_t totalTime = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        times[i] = (totalTime);
        of << times[i] << "\t";
    }
    of << std::endl;
    sort(times.begin(), times.end());
    //erase highest and lowest time
    times.pop_back();
    times.erase(times.begin());


    //calculate mean
    double mean = std::reduce(times.begin(), times.end())/static_cast<double>(times.size());

    //calculate standard deviation
    double standardDev = stdev(times);
    of << mean << "\t"<<standardDev << std::endl;
}
/**
 * Tests bw, fw, bw/r, fw/r, with option for partial and complete
 * 
 */
void CompleteTest(std::ofstream &of, Index &main_index, const std::vector<std::string> &traces, const unsigned char * FMIndex_str, int maxK){


    DoTests(false,false,true, maxK, main_index, traces, FMIndex_str, of);
    of << "\n************************************ **************\n"<<std::endl;
    DoTests(false,true,true, maxK, main_index, traces, FMIndex_str, of);
    of << "\n************************************ **************\n"<<std::endl;
    DoTests(true,false,true, maxK, main_index, traces, FMIndex_str, of);
    of << "\n************************************ **************\n"<<std::endl;
    DoTests(true,true,true, maxK, main_index, traces, FMIndex_str, of);
    of << "\n************************************ **************\n"<<std::endl;
    // Close the file
    of.close();

}
void CompleteTestL2L(std::ofstream &of,  const std::vector<std::string> &traces, int maxK){


    DoTestsL2L(false,false,true, maxK,  traces, of);
    of << "\n************************************ **************\n"<<std::endl;
    DoTestsL2L(false,true,true, maxK,  traces, of);
    of << "\n************************************ **************\n"<<std::endl;
    DoTestsL2L(true,false,true, maxK,  traces, of);
    of << "\n************************************ **************\n"<<std::endl;
    DoTestsL2L(true,true,true, maxK, traces, of);
    of << "\n************************************ **************\n"<<std::endl;
    // Close the file
    of.close();

}
/**
 * Get trace length, optimal value and number of optimal values 
 */
void getOptimalStats(std::ofstream &of, Index &main_index, const std::vector<std::string> &traces, const unsigned char * FMIndex_str, int maxK, bool replace = false){
    ED ed = ED(&main_index, FMIndex_str);
    of << "T No.\tTrace Length \t Optimal Value \t No. of Opt" << std::endl;
    std::vector<int> foundPerK(maxK + 2);

    for(size_t i = 0; i < traces.size(); i++){
        int pK = ed.CheckConf(traces[i], maxK, true, true, replace);
        foundPerK[pK != -1 ? pK : maxK + 1]++; 
        of << i + 1 << "\t";
        of << traces[i].size() << "\t";
        if (pK != -1){
            of << pK << "\t";
            of << ed.getNoOpt(traces[i]) << std::endl;
            // of << traces[i] << std::endl;
            // ed.getAllBest(of);
        }
        
        // of << std::endl;
        
    }
    for (int i = 0; i <= maxK; i++){
        of << i <<":"<< foundPerK[i] << std::endl;
    }
    of<< "NF:" <<foundPerK[maxK + 1] <<std::endl;

}
int main(int argc, char *argv[]){
    int requiredArguments = 5;
    if (argc < requiredArguments){
        std::cerr << "Insufficient number of arguments" << std::endl;
        //show usage
        return EXIT_FAILURE;
    }
    std::cout<< argv[3] << std::endl;
    std::cout<< argv[4] << std::endl;
    std::cout<< argv[5] << std::endl;
    // return 0;
    std::string modelFileName = argv[3];
    std::string logFileName = argv[4];
    std::string resultsFilename = argv[5];
    int K;
    sscanf(argv[6],"%d",&K);
    std::cout<< K << std::endl;

    ConfCheckFMI ccf;
    std::string index_str;
    std::vector<std::string> log;
    try{
        index_str = ccf.load_model(modelFileName);
    }catch (const std::exception& e){
        throw std::runtime_error(e.what());
    }
    try{
        log = ccf.load_log(logFileName);
    }catch (const std::exception& e){
        throw std::runtime_error(e.what());
    }
    std::ofstream file(resultsFilename);
    std::cout << "Index Size:" << index_str.size() << std::endl;
    if (!file.is_open()) {
        std::cerr << "Error opening the file!" << std::endl;
        return 0; 
    }

    // TestIndexes(index_str, file);
    
    const unsigned char * FMIndex_str= reinterpret_cast <const unsigned char *> (index_str.c_str());
    Index idx = makeIndex(FMIndex_str);
    // CompleteTest(file,idx, log, FMIndex_str, K);
    getOptimalStats(file, idx, log,FMIndex_str, K);
    // Log2LogOptimalStats(log,K,file);
    // CompleteTestL2L(file,log,K);
    file.close();
}