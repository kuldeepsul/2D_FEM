#pragma once

#include <iostream>
#include <fstream>
#include <map>
#include <functional>
#include "Model.h"

// Data type
enum class data_type {Node,Frame,BC,Force,Material};


class Input_Reader
{
public:
	data_type current_data{data_type::Node};
	std::map <data_type, std::function<void(model*,std::string&)>> data_handler;
	
	std::ifstream file;

	Input_Reader();


	void read_Node_data(model* p_model,std::string &line);
	void read_Frame_data(model* p_model, std::string& line);
	void read_BC_data(model* p_model, std::string& line);
	void read_Force_data(model* p_model, std::string& line);
	void read_Material_data(model* p_model, std::string& line);
	// Read file
	void read_file(model* p_model,std::string filename);


	// Create Nodes and Elements 
	//model create_model();
};