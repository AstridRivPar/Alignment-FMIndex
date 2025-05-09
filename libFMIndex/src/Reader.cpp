#include <Reader.hpp>

Poco::Logger &Reader::log = Poco::Logger::get("Reader");

Reader::Reader(){}
		
Reader::Reader(std::string model_, std::vector<std::string> queries_){
	model = model_;
	queries = queries_;
}

void Reader::load_model(std::string file_name){
	
	std::string tp;
	std::stringstream ss;
	std::fstream m{file_name};
	int tcounter = 0;
	if(!m.is_open() ){
		std::perror(file_name.c_str());
		exit(-1);
	}
	ss << ",";
	while(getline(m, tp)){
		ss << tp << ",";
		tcounter++;
	}
		
	//ss << tp << "$";
	m.close();
    model = ss.str();
	// std::cout << model << std::endl;
	log.trace("Model Loaded", __FILE__, __LINE__);
	
}

void Reader::load_log_map(std::string file_name){

	std::string line;
	std::fstream p{file_name};
	int tcounter = 0;
	queries.reserve(100);

	if(!p.is_open() ){
		std::perror(file_name.c_str());
		exit(-1);
	}
	while(getline(p, line)){
		std::stringstream linestream(line);
		std::string trace, str_count;
		getline(linestream, trace, ',');
		getline(linestream, str_count);
		int count = stoi(str_count);
		counter.push_back(count);
		queries.push_back(trace);
		tcounter++;

	}
	p.close();
	log.trace("Log Loaded", __FILE__, __LINE__);
}
