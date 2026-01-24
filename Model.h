#pragma once
#include <iostream>
#include <vector>
#include "Material.h"
#include "Matrix.h"

class dof
{
public:
	bool value_known;
	bool is_active{false};
	double value;
	double force;
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
	

	model(int modelid_param) : modelid(modelid_param) {};

	// function to create nodes inside the model
	void create_node(int node_id , double x , double y,double theta);
	void create_frame_element(int element_id, int node_a_id, int node_b_id, Material* mat, double A_param, double I_param);

	Matrix assemble_global_stiffness_matrix();
};