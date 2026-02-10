#pragma once
#include <iostream>
#include <vector>
#include <stack>
#include "Material.h"
#include "Matrix.h"

class dof
{
public:
	bool value_known{false};
	bool is_active{false};
	double value{0};
	double force{0};
	int nodal_dof_id{ 0 };
	const int internal_dof_id;

	dof(int dofid_param) : internal_dof_id(dofid_param)
	{
		// default values 
		this->force = 0;
		this->value = 0;
		this->is_active = true;
		this->value_known = false;
	};
};

class node
{
public:
	std::vector <dof* > nodal_doflist;
	const int nodeid;
    int internal_node_id{0};
	// Data
	const double pos_x;
	const double pos_y;
	const double angle;

	node(int node_id_param, double x, double y, double theta) : nodeid(node_id_param), pos_x(x), pos_y(y), angle(theta) {};

};

class element
{
public:
	int internal_element_id{ 0 };
	std::vector <node* > elemental_nodelist;
	const int elementid;
	Material* mat;
	std::vector <int> dofmap;

	element(int elementid_param) : elementid(elementid_param) {};
	virtual Matrix get_local_stiffness_matrix() = 0;
	virtual void generate_dof_map()= 0;

};

class frame : public element
{
private:
	const int n_nodes{ 3 };
public:
	// Data
	double A{0};
	double I{0};
	Matrix localk;
	node* node_a;
	node* node_b;


	frame(int elementid_param) : element(elementid_param) , localk(n_nodes * 2, n_nodes * 2) {};
	
	Matrix get_local_stiffness_matrix();
	void generate_dof_map();
	
};

class model
{
	const int modelid;
public:
	std::vector <node* > nodelist;
	std::vector <element* > elementlist;
	std::vector <dof* > doflist;
	std::vector <Material*> material_list;
	

	model(int modelid_param) : modelid(modelid_param) {};

	// function to create nodes inside the model
	void create_node(int node_id , double x , double y,double theta);
	void create_frame_element(int element_id, int node_a_id, int node_b_id, int mat_id_param, double A_param, double I_param);
	void create_material(int mat_id_param, double E_param, double Nu_param);

	Matrix assemble_global_stiffness_matrix();
	void prescribe_force(int node_id_param , int dof_id_param , double val);
	void prescribe_dof(int node_id_param, int dof_id_param, double val);

	Matrix get_force_vector();
	//Matrix get_disp_vector();
	Matrix get_reduced_system(Matrix global_k, std::vector <int> &known_dofs);
	Matrix Solve_model(Matrix& stiffness_mat , Matrix& f);
	Matrix generate_full_solution(Matrix& c, std::vector <int> known_dofs);
	//Matrix generate_rection_forces(Matrix& f);

	void export_results(Matrix& sol, std::string outfile_name);
};