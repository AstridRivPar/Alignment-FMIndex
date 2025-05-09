#pragma once
#ifndef READER_H
#define READER_H

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <iostream>
#include <Poco/Logger.h>

class Reader{
	static Poco::Logger &log;

	private:

		std::string model;
		std::vector<std::string> queries;
		std::vector<int> counter;
		// std::vector <int> & variants; //not always used

	public:
		Reader();

		Reader(std::string model_, std::vector<std::string> patterns_);

        /**
         * Reads the file containing the model with the runs as strings and sets
         * the variable model
         * @param file_name the file containing the model
         * 
         */
		void load_model(std::string file_name);
		
		/**
		 * Reads the log where each trace belongs to a variant
		 * @param file_name the file containing the log
		 * @param  vec_var stores the number of variant
		 */
        std::vector<std::string> load_log_var(std::string file_name, std::vector<int> &vec_var);
		/**
		 * Reads the log with the number of traces per unique trace
		 * @param file_name the file containing the log
		 * @param  counter
		 */
		
        void load_log_map(std::string file_name);

		const std::string& get_model() const {return model;}
		const std::vector<std::string>& get_log() const {return queries;}
		const std::vector<int>& get_no_traces() const {return counter;}
};
#endif