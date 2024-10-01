
#include "indexes.hpp"
#include <sdsl/testutils.hpp>
#include "ED.hpp"
#include <ConfCheckFMI.hpp>
#include <algorithm>
#include <chrono>
#include <boost/program_options.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
namespace po = boost::program_options;
namespace logging = boost::log;


void initLogging() {
  auto filter = logging::trivial::severity <= logging::trivial::debug;
  logging::core::get()->set_filter(filter);
}

index_bidirectional_waveletindex<> BuildIndex(const unsigned char * FMIndex_str){
    
    index_bidirectional_waveletindex<> main_index = index_bidirectional_waveletindex(FMIndex_str);
    BOOST_LOG_TRIVIAL(trace) << "SUCCESS BUILDING INDEX";
    // return FMIndex_str;
    return main_index;
}

void printAll(std::ofstream& of, const std::map<std::string, std::vector<Result>> &results){
    for (auto const & result: results){
        if (!result.second.empty()){
            of <<"Trace: \t" << result.first << "\tk: "<< result.second[0].k << "\tNo. opt: "<< result.second.size()<<std::endl;
            for (const auto &element: result.second){
                of << *element.match << "\t";
                for (auto c = element.alignment->begin(); c!=element.alignment->end(); c++){
                    of << *c;
                }  
                of << std::endl;
            }
        }else{
            of <<"Trace: \t" << result.first << "No matches" <<std::endl;
        }
    }
}
int main( int argc, char *argv[]){
    initLogging();
    po::options_description desc("\nMandatory arguments marked with '*'.\n"
    "Invocation: <program> --backward <optional> -- partial <optional> <index> <queries> <k>\n Allowed options");
    desc.add_options()
        ("b", "do the search backwards")
        ("p", "match partial traces")
        ("index", po::value<std::string>()->required(), "* Path to file to build index")
        ("queries", po::value<std::string>()->required(), "* Path to file with queries")
        ("k", po::value<int>()->required(), "* Maximum number of mismatches")
        ("output", po::value<std::string>()->required(), "* Path to outputfile")
        ;
    po::positional_options_description pos_desc;
    pos_desc.add("index", 1);
    pos_desc.add("queries", 1);
    pos_desc.add("k", 1);
    pos_desc.add("output", 1); //*********** 
    po::variables_map vm;
    
    try{
        po::store(po::command_line_parser(argc, argv).options(desc).positional(pos_desc).run(), vm);
        po::notify(vm);
       
    }
    catch(po::error& e){
        std::cout << "Error: " << e.what() << "\n";
        std::cout << desc << "\n";
        exit(EXIT_FAILURE);

    }
    int k = 0;
    bool fw = true, complete = true, replace = false;
    std::string queries_file, index_file, out_filename;

    
    k = vm["k"].as<int>();
    index_file = vm["index"].as<std::string>();
    queries_file = vm["queries"].as<std::string>();
    if (vm.count("b")) fw = false;
    if (vm.count("p")) complete = false;
    std::ofstream outFile;
    if (vm.count("output")) {
        out_filename = vm["output"].as<std::string>();
        outFile.open(out_filename);
        if(!outFile.is_open()){
            std::perror(out_filename.c_str());
            exit(EXIT_FAILURE); 
        }
    }
    BOOST_LOG_TRIVIAL(trace) << "\nParameters: \nk: " << k <<"\nIndex file: "<< index_file <<
                                "\nQueries file: " << queries_file <<
                                "\nForward: " << fw << "\nComplete: " << complete;
    
    
    ConfCheckFMI ccf;
    
    auto index_str = ccf.load_model(index_file);
    auto queries = ccf.load_log(queries_file);
   
    BOOST_LOG_TRIVIAL(trace) << "Index Size:" << index_str.size();
    const unsigned char * FMIndex_str= reinterpret_cast <const unsigned char *> (index_str.c_str());
    auto FMIndex = BuildIndex(FMIndex_str);
    
     

    std::chrono::high_resolution_clock::time_point start, end;
    start = std::chrono::high_resolution_clock::now(); 
    auto as = ED(&FMIndex, FMIndex_str);
    auto results = as.CompleteAligment(queries,k, fw, complete,replace);
    end = std::chrono::high_resolution_clock::now(); 
    
    int totalTime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
    
    printAll(outFile, results);
    for (int i = 0; i <= k; i++){
        std::cout<< i <<": "<< as.getFoundPerK()[i] << std::endl;
    }
    std::cout<< "Not found: "<< as.getFoundPerK()[k + 1] << std::endl;
    std::cout<< "Total number of patterns: "<< queries.size() << std::endl;
    std::cout<< "Total time (ms): "<< totalTime << std::endl;
    std::cout<< "Average time per pattern (ms): "<< totalTime/queries.size()<< std::endl;
    std::cout<< "Total time (s): "<< totalTime/(1000.0)<< std::endl;
    outFile.close();
}