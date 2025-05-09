
#include <sstream>
#include <iostream>
#include <Poco/Util/Option.h>
#include <Poco/Util/OptionSet.h>
#include <Poco/Util/Application.h>
#include <Poco/Util/HelpFormatter.h>
#include <Poco/Util/OptionCallback.h>
// #include <vector>
#include "Bi_FM_Index.hpp"
#include "Reader.hpp"
#include <algorithm>
#include "ED.hpp"
// #include "config.h"

#include <string.h>

using Poco::Util::Application;
using Poco::Util::Option;
using Poco::Util::OptionSet;
using Poco::Util::HelpFormatter;
using Poco::Util::OptionCallback;

class FMIndex_App: public Application{

	Poco::Logger &log = Poco::Logger::get("FMIndex_app");

	private:

		bool helpRequested;
		std::string direction; 				//default true (backward)
		std::string complete;
		std::string index_file;
        std::string queries_file;
		std::string output_file;
		bool d = true;
        int k;
		bool c = true;

	public:
        FMIndex_App(): helpRequested{false}{}

	protected:	

		void initialize(Application& self) {
			loadConfiguration();
			Application::initialize(self);
		}

		void uninitialize() {
			Application::uninitialize();
		}

		void reinitialize(Application& self) {
			Application::reinitialize(self);
		}

        void defineOptions(OptionSet& options){
            Application::defineOptions(options);
        

			options.addOption(
					Option("help", "h", "Display this help information")
					.callback(OptionCallback<FMIndex_App>(this, &FMIndex_App::handleHelp)));

			options.addOption(
					Option("forward", "f", "forward matching")
					.callback(OptionCallback<FMIndex_App>(this, &FMIndex_App::set_direction)));

            options.addOption(
					Option("backward", "b", "backward matching(default)")
					.callback(OptionCallback<FMIndex_App>(this, &FMIndex_App::set_direction)));

			options.addOption(
						Option("complete", "c", "matches complete traces(default)")
						.callback(OptionCallback<FMIndex_App>(this, &FMIndex_App::set_complete)));
			
			options.addOption(
				Option("partial", "p", "matches partial traces")
				.callback(OptionCallback<FMIndex_App>(this, &FMIndex_App::set_complete)));

			options.addOption(
					Option("k", "k", "maximum number of mismatches")
                    .required(true)
                    .argument("<int>")
					.callback(OptionCallback<FMIndex_App>(this, &FMIndex_App::set_maxK)));

			options.addOption(
					Option("index_file", "i", "input file name to build index")
					.required(true)
					.argument("<file_name>")
					.callback(OptionCallback<FMIndex_App>(this, &FMIndex_App::set_index_file)));
            
            options.addOption(
					Option("query_file", "q", "file with the queries")
					.required(true)
					.argument("<file_name>")
					.callback(OptionCallback<FMIndex_App>(this, &FMIndex_App::set_query_file)));

			options.addOption(
					Option("output_file", "o", "output file name")
                    .required(true)
					.argument("<file_name>")
					.callback(OptionCallback<FMIndex_App>(this, &FMIndex_App::set_output_file)));


		}

		void set_complete(const std::string& name, const std::string& value){
			complete = name;
			c = complete == "complete"? true: false;
		}
		void set_index_file(const std::string& name, const std::string& value) {
			index_file = value;
		}

		void set_output_file(const std::string& name, const std::string& value) {
			output_file = value;
		}

        void set_query_file(const std::string& name, const std::string& value) {
			queries_file = value;
		}

		void set_direction(const std::string& name, const std::string& value) { 
			direction = name;
			d = direction == "forward"? true:false;
		}

		void set_maxK(const std::string& name, const std::string& value) { 
			k = stoi(value);
		}

		void handleHelp(const std::string& name, const std::string& value) {
			helpRequested = true;
			displayHelp();
			stopOptionsProcessing();
		}

		void displayHelp() {
			HelpFormatter helpFormatter(options());
			helpFormatter.setCommand(commandName());
			helpFormatter.setUsage("<options>");
			helpFormatter.setHeader("A command line interface (cli) application for matchin runs with traces from a log.");
			helpFormatter.format(std::cout);
		}

        int main(const std::vector<std::string> &args){
            if (!helpRequested){

				log.information(std::string{"***** Running alignments "});

                log.information("-----Application options-----" , __FILE__,__LINE__);
                log.information("Direction  : " + direction, __FILE__,__LINE__);
				log.information("K          : " + std::to_string(k), __FILE__,__LINE__);
				log.information("Is complete:" +  complete, __FILE__,__LINE__);
                log.information("Index file : " + index_file, __FILE__,__LINE__);
                log.information("Query file : " + queries_file, __FILE__,__LINE__);
                log.information("Output file: " + output_file, __FILE__,__LINE__);
                
				
				std::ofstream outFile (output_file);
				if(!outFile.is_open()){
					std::perror(output_file.c_str());
					exit(EXIT_FAILURE); 
				}
				
				//Load model and log
				
				Reader ccf = Reader();
				ccf.load_model(index_file);
				ccf.load_log_map(queries_file);

				// Build index
				log.trace("Building index", __FILE__,__LINE__);
				Bi_FM_Index index(ccf.get_model());

				//Match
				ED ed(&index, &ccf.get_model());
				ed.CompleteAligment(ccf.get_log(), k, d, c);
				// for (int i = 0; i <= k; i++){
				// 	std::cout<< i <<": "<< ed.getFoundPerK()[i] << std::endl;
				// }

				ed.printAll(outFile, ccf);


			}
            return Application::EXIT_OK;
        }

};

POCO_APP_MAIN(FMIndex_App)
