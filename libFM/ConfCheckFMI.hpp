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
		std::pair<std::string, std::vector<std::pair<size_t, size_t>>> load_model_var(std::string file_name);
		std::vector<std::string> load_log(std::string file_name);
		std::vector<std::string> load_log(std::string file_name, std::unordered_map<std::string, size_t> &traces_with_id);
		std::vector<std::string> load_log_var(std::string file_name, std::unordered_map<std::string, std::vector<uint8_t>> &vec_var);
		
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
std::pair<std::string, std::vector<std::pair<size_t, size_t>>> ConfCheckFMI::load_model_var(std::string file_name){
	std::string tp;
	std::stringstream ss;
	std::fstream model{file_name};
	std::vector<std::pair<size_t, size_t>> var_index;
	var_index.reserve(5);

	if(!model.is_open() ){
		std::perror(file_name.c_str());
		exit(-1);
	}
	
	ss << ",";
	size_t counter_size = 1;
	size_t st_idx = 0, end_idx = 0;
	while(getline(model, tp)){
		if (tp == "#"){
			end_idx = counter_size;
			var_index.emplace_back(st_idx, end_idx);
			st_idx =  end_idx;
			ss << ",";
			counter_size ++;
			continue;
		}
		ss << tp << ",";
		counter_size+=tp.size()+1;
		
	}
		
	//ss << tp << "$";
	model.close();
	//std::cout << ss.str() << std::endl;
	return std::make_pair(ss.str(), var_index);

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
/**
 * Reads with groud truth to which variant it belongs to
 */
std::vector<std::string> ConfCheckFMI::load_log_var(std::string file_name, std::unordered_map<std::string, std::vector<uint8_t>> &vec_var){

	std::string line;
	std::fstream patterns{file_name};
	std::vector<std::string> patternsv;
	patternsv.reserve(100);

	if(!patterns.is_open() ){
		std::perror(file_name.c_str());
		exit(-1);
	}
	while(getline(patterns, line)){
		
	// 	results[query] = vector<Result>();
    // sols = &results[query];
		std::stringstream linestream(line);
		std::string trace, str_variant;
		getline(linestream, trace, ',');
		getline(linestream, str_variant, ',');
		uint8_t variant = stoi(str_variant);

		std::vector<uint8_t> *variants;
		if (vec_var.find(trace) == vec_var.end()){

			vec_var[trace] = vector<uint8_t>();
			variants = &vec_var[trace];
			variants->reserve(5);
			// patternsv.push_back(tp);
		}
		else{
			variants = &vec_var[trace];
		}
		
		variants->push_back(variant);
		patternsv.push_back(trace);
	}
	patterns.close();

	return patternsv;

}
