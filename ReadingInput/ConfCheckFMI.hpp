#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <unordered_map>
class ConfCheckFMI{

	private:

		std::string model;
		std::vector<std::string> patterns;

	public:
		ConfCheckFMI();
		ConfCheckFMI(std::string model_, std::vector<std::string> patterns_);

		std::string load_model(std::string file_name);
		std::vector<std::string> load_log(std::string file_name);
		std::vector<std::string> load_log(std::string file_name, std::unordered_map<std::string, size_t> &traces_with_id);
};

ConfCheckFMI::ConfCheckFMI(){}
		
ConfCheckFMI::ConfCheckFMI(std::string model_, std::vector<std::string> patterns_){
	model = model_;
	patterns = patterns_;
}

std::string ConfCheckFMI::load_model(std::string file_name){
	std::string tp;
	std::stringstream ss;
	std::fstream model{file_name};

	if(!model.is_open() ){
		std::perror(file_name.c_str());
		exit(-1);
	}
	ss << ",";
	while(getline(model, tp))
		ss << tp << ",";
	//ss << tp << "$";
	model.close();
	//std::cout << ss.str() << std::endl;
	return ss.str();

}

std::vector<std::string> ConfCheckFMI::load_log(std::string file_name){

	std::string tp;
	std::fstream patterns{file_name};
	std::vector<std::string> patternsv;

	if(!patterns.is_open() ){
		std::perror(file_name.c_str());
		exit(-1);
	}

	while(getline(patterns, tp))
		patternsv.push_back(tp);
	patterns.close();

	return patternsv;

}

std::vector<std::string> ConfCheckFMI::load_log(std::string file_name, std::unordered_map<std::string, size_t> &traces_with_id){

	std::string tp;
	std::fstream patterns{file_name};
	std::vector<std::string> patternsv;

	if(!patterns.is_open() ){
		std::perror(file_name.c_str());
		exit(-1);
	}
	size_t counter = 0;
	while(getline(patterns, tp)){
		patternsv.push_back(tp);
		traces_with_id[tp] = counter++;
	}
	patterns.close();

	return patternsv;

}
