#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include "Matrix.h"
#include "Model.h"
#include "Material.h"

int main()
{
	/*
	Matrix Trial(5,5);
	Trial.data = {2,4,5,1,2, 4,1,5,2,7, 8,4,7,3,0, 0,4,7,3,1, 5,6,1,7,3};
	Matrix f(5, 1);
	f.data = {5,1,10,2,6};

	Trial.print_matrix();
	*/

	
	Material mat_01;
	mat_01.E = 10e8;
	mat_01.nu = 0.28;
	model model_01(0);
	model_01.create_node(1, 0, 0, 0);
	model_01.create_node(2, 10, 0, 0);
	model_01.create_node(3, 20, 0, 0);

	model_01.create_frame_element(1, 1, 2, &mat_01, 0.1, 0.0001);
	model_01.create_frame_element(2, 2, 3, &mat_01, 0.1, 0.0001);

	Matrix global_k = model_01.assemble_global_stiffness_matrix();
	model_01.prescribe_dof(1, 0, 0);
	model_01.prescribe_dof(1, 1, 0);
	model_01.prescribe_dof(1, 2, 0);
	model_01.prescribe_force(3,2,1000);

	global_k.print_matrix();

	std::stack <int> known_dofs;
	Matrix reduced_global_k = model_01.get_reduced_system(global_k,known_dofs);
	std::cout << "Reduced system" << std::endl;
	reduced_global_k.print_matrix();
	Matrix f = model_01.get_force_vector();
	std::cout << "___________________________" << std::endl;
	f.print_matrix();
	std::cout << "___________________________" << std::endl;
	Matrix c = model_01.Solve_model(reduced_global_k,f);
	c.print_matrix();
	Matrix x = model_01.generate_full_solution(c);



	/*
	int rows = Trial.rows;
	int cols = Trial.cols;
	Matrix a(rows, 1);
	Matrix b(rows, 1);
	Matrix c(rows, 1);
	

	Matrix P(rows,cols), L(rows,cols), U(rows,cols);
	Trial.PLUDecomposition(P,L,U);
	std::cout << "PLU decomposition complete ! " << std::endl;
	P.print_matrix();
	std::cout << "-----------------------------------------" << std::endl;
	L.print_matrix();
	std::cout << "-----------------------------------------" << std::endl;
	U.print_matrix();
	std::cout << "-----------------------------------------" << std::endl;
	f.print_matrix(); 
	std::cout << "--------------------------------------------" << std::endl;
	std::cout << "Solving P" << std::endl;
	a = Trial.Solve_linear_system(P, f,true);
	a.print_matrix();
	std::cout << "Solving L" << std::endl;
	b = Trial.Solve_linear_system(L, a, true);
	b.print_matrix();
	std::cout << "Solving U" << std::endl;
	c = Trial.Solve_linear_system(U, b, false);
	c.print_matrix();
	*/

}