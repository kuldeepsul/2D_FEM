#pragma once
#include "Model.h"
#include <cmath>

Matrix frame::get_local_stiffness_matrix()
{
	// generating the values required for calculation
	double lx = this->elemental_nodelist[1]->pos_x - this->elemental_nodelist[0]->pos_x;
	double ly = this->elemental_nodelist[1]->pos_y - this->elemental_nodelist[0]->pos_y;
	double l = sqrt((lx * lx) + (ly * ly));
	double c = lx / l;
	double s = ly / l;
	double ka = (this->A * this->mat->E) / l;
	double kb = (this->I * this->mat->E) / (l*l*l);

	// generating required matrices
	Matrix T(this->localk.rows , this->localk.cols);
	Matrix elementk(this->localk.rows, this->localk.cols);

	T.data = { c , -s , 0 , 0 , 0 , 0 ,
			   s ,  c , 0 , 0 , 0 , 0 ,
			   0 ,  0 , 1 , 0 , 0 , 0 ,
			   0 ,  0 , 0 , c ,-s , 0 ,
			   0 ,  0 , 0 , s , c , 0 ,
			   0 ,  0 , 0 , 0 , 0 , 1 };

	elementk.data = { ka ,  0          , 0				, -ka, 0		   , 0				,
					   0 ,  12 * kb    , 6 * l * kb     , 0  , -12 * kb    , 6 * l * kb     ,
					   0 ,  6 * l * kb , 4 * l * l * kb , 0  , -6 * l * kb , 2 * l * l * kb ,
					 -ka ,  0          , 0              , ka , 0		   , 0				,
					   0 ,  -12 * kb   , -6 * l * kb 	, 0  , 12 * kb     , -6 * l * kb    ,
					   0 , 6 * l * kb  , 2 * l * l * kb , 0  ,-6 * l * kb  , 4 * l * l * kb };
	this->localk = T * elementk;
	this->localk = T.Transpose() * localk;
	return localk;
};

void model::create_node(int node_id, double x, double y, double theta)
{
	// genrating node and pushing into the global node list
	node* n = new node(node_id, x, y, theta);
	n->internal_node_id = size(this->nodelist);
	this->nodelist.push_back(n);

	// generating three dof for the node and pushing into global dof list
	dof* p_dof_01 = new dof(size(this->doflist)); 
	this->doflist.push_back(p_dof_01);
	dof* p_dof_02 = new dof(size(this->doflist));
	this->doflist.push_back(p_dof_02);
	dof* p_dof_03 = new dof(size(this->doflist));
	this->doflist.push_back(p_dof_03);

	// pushing dof into nodal dof list
	p_dof_01->nodal_dof_id = size(n->nodal_doflist);
	n->nodal_doflist.push_back(p_dof_01);
	p_dof_02->nodal_dof_id = size(n->nodal_doflist);
	n->nodal_doflist.push_back(p_dof_02);
	p_dof_03->nodal_dof_id = size(n->nodal_doflist);
	n->nodal_doflist.push_back(p_dof_03);

	
}

void model::create_frame_element(int element_id, int node_a_id, int node_b_id, Material* mat_param, double A_param, double I_param)
{
	frame* new_element = new frame(element_id);
	new_element->mat = mat_param;
	new_element->A = A_param;
	new_element->I = I_param;

	bool node_a_found{false};
	bool node_b_found{false};

	int initial_node_count = size(this->nodelist);

	for (int i{ 0 }; i < initial_node_count; i++)
	{
		if (!node_a_found)
		{
			if (this->nodelist[i]->nodeid == node_a_id)
			{
				new_element->elemental_nodelist.push_back(nodelist[i]);
				new_element->node_a = nodelist[i];
				node_a_found = true;
			}
			else if (i + 1 == size(this->nodelist) && (!node_a_found))
			{
				std::cout << " Error : Node Id not found : " << node_a_id << std::endl;
			}
		}

		if (!node_b_found)
		{
			if (this->nodelist[i]->nodeid == node_b_id)
			{
				new_element->elemental_nodelist.push_back(nodelist[i]);
				new_element->node_b = nodelist[i];
				node_b_found = true;
			}
			else if (i + 1 == size(this->nodelist) && (!node_b_found))
			{
				std::cout << " Error : Node Id not found : " << node_a_id << std::endl;
			}
		}

	}
	new_element->get_local_stiffness_matrix();
	this->elementlist.push_back(new_element);


}