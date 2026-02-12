#include "Input_Reader.h"
#include <filesystem>

Input_Reader::Input_Reader()
{
	current_data = data_type::Node;
	data_handler[data_type::Node] = [this](model* m,std::string &line) { read_Node_data(m,line); };
	data_handler[data_type::Frame] = [this](model* m, std::string& line) { read_Frame_data(m,line); };
	//data_handler[data_type::Quad] = [this](model* m, std::string& line) { read_Quad_data(m, line); };
	data_handler[data_type::BC] = [this](model* m, std::string& line) { read_BC_data(m,line); };
	data_handler[data_type::Force] = [this](model* m, std::string& line) { read_Force_data(m,line); };
	data_handler[data_type::Material] = [this](model* m, std::string& line) { read_Material_data(m, line); };
};

void Input_Reader::read_Node_data(model* p_model,std::string &line)
{
	
	std::stringstream s(line);
	std::string value;

	double node_id, x, y;
	std::getline(s, value, ',');
	node_id = std::stoi (value);
	std::getline(s, value, ',');
	x = std::stod(value);
	std::getline(s, value, ',');
	y = std::stod(value);

	// For debugging
	std::cout << "Node_ID: " << node_id << " X:" << x << " Y:" << y << std::endl;

	// Creating Instances
	p_model->create_node(node_id,x,y,0.0);

}
void Input_Reader::read_Frame_data(model* p_model, std::string& line)
{
	
	std::string value;
	std::stringstream s(line);

	double element_id, n1, n2, A, I, mat_id;
	std::getline(s, value, ',');
	element_id = std::stoi(value);
	std::getline(s, value, ',');
	n1 = std::stoi(value);
	std::getline(s, value, ',');
	n2 = std::stoi(value);
	std::getline(s, value, ',');
	A = std::stod(value);
	std::getline(s, value, ',');
	I = std::stod(value);
	std::getline(s, value, ',');
	mat_id = std::stoi(value);

	// For debugging
	std::cout << "Element_ID: " << element_id << " N1:" << n1 << " N2:" << n2  << " A:" << A << " I:" << I << " Materail_ID:" << mat_id << std::endl;

	// Creating Instance
	p_model->create_frame_element(element_id,n1,n2,mat_id,A,I);
}
void Input_Reader::read_BC_data(model* p_model, std::string& line)
{
	
	std::string value;
	std::stringstream s(line);

	double node_id, dof_id, val;
	std::getline(s, value, ',');
	node_id = std::stoi(value);
	std::getline(s, value, ',');
	dof_id = std::stoi(value);
	std::getline(s, value, ',');
	val = std::stod(value);

	// For debugging
	std::cout << "Node_ID: " << node_id << " DOF:" << dof_id << " Value:" << val << std::endl;

	// Prescribing Boundary Condition
	p_model->prescribe_dof(node_id,dof_id,val);


}
void Input_Reader::read_Force_data(model* p_model, std::string& line)
{
	
	std::string value;
	std::stringstream s(line);

	double node_id, dof_id, val;
	std::getline(s, value, ',');
	node_id = std::stoi(value);
	std::getline(s, value, ',');
	dof_id = std::stoi(value);
	std::getline(s, value, ',');
	val = std::stod(value);

	// For debugging
	std::cout << "Node_ID: " << node_id << " DOF:" << dof_id << " Magnitude:" << val << std::endl;

	// Prescribing Force
	p_model->prescribe_force(node_id,dof_id,val);

}

void Input_Reader::read_Material_data(model* p_model, std::string& line)
{

	std::string value;
	std::stringstream s(line);

	double mat_id, E, nu;
	std::getline(s, value, ',');
	mat_id = std::stoi(value);
	std::getline(s, value, ',');
	E = std::stod(value);
	std::getline(s, value, ',');
	nu = std::stod(value);

	// For debugging
	std::cout << "Material_ID: " << mat_id << " Young's Modulus:" << E << " Poission's Ratio:" << nu << std::endl;

	// Creating Instances
	p_model->create_material(mat_id,E,nu);
}

void Input_Reader::read_file(model* p_model,std::string filename)
{
	
	this->file.open(filename);
	if (file.is_open())
	{
		std::cout << "Input file opened sucessfully" << std::endl;
	}
	std::string line;

	while (std::getline(file, line))
	{
		if (line == "*Node")
		{
			this->current_data = data_type::Node;
			std::cout << "Node ----------------------------------------------" << std::endl;
		}
		else if (line == "*Frame")
		{
			this->current_data = data_type::Frame;
			std::cout << "Frame Element -------------------------------------" << std::endl;
		}
		else if (line == "*BC")
		{
			this->current_data = data_type::BC;
			std::cout << "Boundary Condition --------------------------------" << std::endl;
		}
		else if (line == "*Force")
		{
			this->current_data = data_type::Force;
			std::cout << "Force ---------------------------------------------" << std::endl;
		}
		else if (line == "*Material")
		{
			this->current_data = data_type::Material;
			std::cout << "Material ------------------------------------------" << std::endl;
		}
		else
		{
			this->data_handler[this->current_data](p_model,line);
		}
	}

	std::cout << "Input Reader Finished -------------------------------------" << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
}