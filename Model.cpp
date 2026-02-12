#pragma once
#include "Model.h"
#include <cmath>
#include <fstream>

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

	T.data = { c ,  s , 0 , 0 , 0 , 0 ,
			  -s ,  c , 0 , 0 , 0 , 0 ,
			   0 ,  0 , 1 , 0 , 0 , 0 ,
			   0 ,  0 , 0 , c , s , 0 ,
			   0 ,  0 , 0 ,-s , c , 0 ,
			   0 ,  0 , 0 , 0 , 0 , 1 };


	elementk.data = { ka ,  0          , 0				, -ka, 0		   , 0				,
					   0 ,  12 * kb    , 6 * l * kb     , 0  , -12 * kb    , 6 * l * kb     ,
					   0 ,  6 * l * kb , 4 * l * l * kb , 0  , -6 * l * kb , 2 * l * l * kb ,
					 -ka ,  0          , 0              , ka , 0		   , 0				,
					   0 ,  -12 * kb   , -6 * l * kb 	, 0  , 12 * kb     , -6 * l * kb    ,
					   0 , 6 * l * kb  , 2 * l * l * kb , 0  ,-6 * l * kb  , 4 * l * l * kb };

	this->localk = elementk * T;

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
	p_dof_01->nodal_dof_id = 1;
	n->nodal_doflist.push_back(p_dof_01);
	p_dof_02->nodal_dof_id = 2;
	n->nodal_doflist.push_back(p_dof_02);
	p_dof_03->nodal_dof_id = 3;
	n->nodal_doflist.push_back(p_dof_03);

	
}


void model::create_frame_element(int element_id, int node_a_id, int node_b_id,int mat_id_param, double A_param, double I_param)
{
	frame* new_element = new frame(element_id);
	new_element->A = A_param;
	new_element->I = I_param;

	bool node_a_found{false};
	bool node_b_found{false};
	bool material_found{false};

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
	////////////////////////////////////////////////////////////////
	// Assigning Material
	for (int i{ 0 }; i < size(this->material_list); i++)
	{
		if (this->material_list[i]->M_id == mat_id_param)
		{
			new_element->mat = material_list[i];
			material_found = true;
		}
		else
		{
			continue;
		}
		
	}
	if (!material_found)
	{
		std::cout << "Error : Material Id not found : " << element_id << std::endl;
	}
	//////////////////////////////////////////////////////////////////
	new_element->get_local_stiffness_matrix();
	new_element->generate_dof_map();
	this->elementlist.push_back(new_element);


}


//////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////			Q4 Functions						//////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

Matrix Quad::get_shape_func_der(double xi, double eta)
{
	Matrix sf_der{2,4};

	sf_der.data = { -(1 - eta) / 4.0	,(1 - eta) / 4.0	,(1 + eta) / 4.0	,-(1 + eta) / 4.0,
					-(1 - xi) / 4.0 	,-(1 + xi) / 4.0	,(1 + xi) / 4.0		,(1 - xi) / 4.0 };

	return sf_der;
}

Matrix Quad::get_jacobian_matrix(double xi, double eta)
{
	Matrix sf_der = this->get_shape_func_der(xi, eta);

	return sf_der * this->corners;
}

Matrix Quad::get_B_matrix(Matrix& sf_der, Matrix& jacobian)
{
	Matrix I_jacobian = jacobian.Invert_Jacobain(jacobian);

	Matrix N1_der_n{ 2,1 };
	Matrix N2_der_n{ 2,1 };
	Matrix N3_der_n{ 2,1 };
	Matrix N4_der_n{ 2,1 };

	N1_der_n.data = { sf_der(0,0) , sf_der(1,0) };
	N2_der_n.data = { sf_der(0,1) , sf_der(1,1) };
	N3_der_n.data = { sf_der(0,2) , sf_der(1,2) };
	N4_der_n.data = { sf_der(0,3) , sf_der(1,3) };

	Matrix N1_der_g = I_jacobian * N1_der_n;
	Matrix N2_der_g = I_jacobian * N2_der_n;
	Matrix N3_der_g = I_jacobian * N3_der_n;
	Matrix N4_der_g = I_jacobian * N4_der_n;

	Matrix B{ 3,8 };
	B.data = { N1_der_g(0,0) ,		0		,N2_der_g(0,0)	,		0		,N3_der_g(0,0) ,		0	  ,N4_der_g(0,0) ,		0
			 ,     0	     ,N1_der_g(1,0) ,		0		,N2_der_g(1,0)	,		0	   ,N3_der_g(0,0) ,		0		 ,N4_der_g(0,0)
			 ,N1_der_g(1,0)  ,N1_der_g(0,0) ,N2_der_g(1,0)  ,N2_der_g(0,0)  ,N3_der_g(1,0) ,N3_der_g(0,0) ,N4_der_g(1,0) ,N4_der_g(0,0)};


}

Matrix Quad::get_local_stiffness_matrix()
{
	Matrix Ke{ 8,8 };
	std::vector <std::vector <double>> samples;
	samples = { {-0.57735,-0.57735},{0.57735,-0.57735},{0.57735,0.57735},{-0.57735,0.57735} };

	for (int i{ 0 }; i < size(samples); i++)
	{
		Matrix sf_der = this->get_shape_func_der(samples[i][0], samples[i][1]);
		Matrix J = this->get_jacobian_matrix(samples[i][0], samples[i][1]);
		Matrix B = this->get_B_matrix(sf_der, J);

		Matrix sample_matrix = this->mat->D_stress * B;
		Matrix sample_matrix = B.Transpose() * localk;

		Ke = Ke + sample_matrix;


	}
}





void model::create_material(int mat_id_param, double E_param, double Nu_param)
{
	Material* mat = new Material(mat_id_param,E_param,Nu_param);
	this->material_list.push_back(mat);
}

Matrix model::assemble_global_stiffness_matrix()
{
	size_t n = this->doflist.size();
	Matrix gk((int(n)),(int(n)));
	double val = 0;

	for (int i{ 0 }; i < size(this->elementlist); i++)
	{
		Matrix temp = elementlist[i]->get_local_stiffness_matrix();
		for (int j{ 0 }; j < size(elementlist[i]->dofmap); j++)
		{
			for (int k{ 0 }; k < size(elementlist[i]->dofmap); k++)
			{
				gk((elementlist[i]->dofmap[j]), (elementlist[i]->dofmap[k])) += temp(j, k);
			}
		}
	}
	return gk;
}

void frame::generate_dof_map()
{
	for (int i{ 0 }; i < size(this->elemental_nodelist); i++)
	{
		for (int k{ 0 }; k < size(this->elemental_nodelist[i]->nodal_doflist); k++)
		{

			this->dofmap.push_back(this->elemental_nodelist[i]->nodal_doflist[k]->internal_dof_id);
			
		}
	}
}

void model::prescribe_force(int node_id_param, int dof_id_param, double val)
{
	for (int i{ 0 }; i < size(this->nodelist); i++)
	{
		if (node_id_param == this->nodelist[i]->nodeid)
		{
			for (int k{ 0 }; k < size(this->nodelist[i]->nodal_doflist); k++)
			{
				if (dof_id_param == this->nodelist[i]->nodal_doflist[k]->nodal_dof_id)
				{
					this->nodelist[i]->nodal_doflist[k]->force = val;
				}
			}
		}
	}

}
void model::prescribe_dof(int node_id_param, int dof_id_param, double val)
{
	for (int i{ 0 }; i < size(this->nodelist); i++)
	{
		if (node_id_param == this->nodelist[i]->nodeid)
		{
			for (int k{ 0 }; k < size(this->nodelist[i]->nodal_doflist); k++)
			{
				if (dof_id_param == this->nodelist[i]->nodal_doflist[k]->nodal_dof_id)
				{
					this->nodelist[i]->nodal_doflist[k]->value = val;
					this->nodelist[i]->nodal_doflist[k]->value_known = true;
				}
			}
		}
	}
}

Matrix model::get_force_vector()
{
	
	Matrix f (0,0);

	for (int i{ 0 }; i < size(this->doflist); i++)
	{
		if (!this->doflist[i]->value_known)
		{
			f.data.push_back(this->doflist[i]->force);
		}
	}
	f.rows = size(f.data);
	f.cols = 1;
	return f;
}

Matrix model::get_reduced_system(Matrix global_k, std::vector <int> &known_dofs)
{
	int val = global_k.cols;
	
	
	for (int i{ 0 }; i < size(this->doflist); i++)
	{
		if (this->doflist[i]->value_known)
		{
			known_dofs.push_back(i);
		}
	}

	std::vector <int> temp = known_dofs;
	int iters = size(temp);
	for (int j{ 0 }; j < iters; j++)
	{
		global_k.remove_row(temp.back());
		temp.pop_back();
	}

	return global_k;
}

Matrix model::generate_full_solution(Matrix& c, std::vector <int> known_dofs)
{
	
	int iter = size(this->doflist);
	Matrix x(iter , 1);

	std::vector <int> temp = known_dofs;
	for (int i{ 0 }; i < iter; i++)
	{
		if (size(temp) != 0)
		{
			if (i == temp.front())
			{
				x(i, 0) = 0;
				temp.erase(temp.begin());
				continue;
			}
			
		}

		x(i, 0) = c.data.front();
		c.data.erase(c.data.begin());
		
	}
	return x;
};
//Matrix generate_rection_forces(Matrix& f);

Matrix model::Solve_model(Matrix& stiffness_mat, Matrix& f)
{
	// Generating reduced stiffness matrix
	std::vector <int> known_dofs;
	Matrix reduced_stiffness_mat = this->get_reduced_system(stiffness_mat,known_dofs);

	std::cout << "Reduced System--------------------------------------------------" << std::endl;
	reduced_stiffness_mat.print_matrix();

	// PLU decompostion of stiffness matrix
	int rows = reduced_stiffness_mat.rows;
	int cols = reduced_stiffness_mat.cols;
	Matrix P(rows,cols);
	Matrix L(rows, cols);
	Matrix U(rows, cols);
	// create intermediate vectors
	Matrix a(rows, 1);
	Matrix b(rows, 1);
	Matrix c(rows, 1);

	reduced_stiffness_mat.PLUDecomposition(P,L,U);

	// [P]{a} = {f}
	a = reduced_stiffness_mat.Solve_linear_system(P, f, true);
	// [L]{b} = {a}
	b = reduced_stiffness_mat.Solve_linear_system(L, a, true);
	// [U]{c} = {b}
	c = reduced_stiffness_mat.Solve_linear_system(U, b , false);

	Matrix x = this->generate_full_solution(c,known_dofs);
	return x;
}

void model::export_results(Matrix& sol, std::string outfile_name)
{
	std::fstream res_file;
	res_file.open(outfile_name,std::ios::out | std::ios::app);
	int count{ 0 };
	
	for (int i{ 0 }; i < size(this->nodelist); i++)
	{
		for (int k{ 0 }; k < size(this->nodelist[i]->nodal_doflist); k++)
		{
			res_file << this->nodelist[i]->nodeid << "," << this->nodelist[i]->nodal_doflist[k]->nodal_dof_id << "," << sol(count, 0) << std::endl;
			count++;
		}
		
	}
}

