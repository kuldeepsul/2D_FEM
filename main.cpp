#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include "Matrix.h"
#include "Model.h"
#include "Material.h"
#include "Input_Reader.h"

int main()
{	/*
	model* p_model = new model(1);
	
	////////////////////////////////////////////////////////////////////////////////////
	// Data Input from .inp file

	
	Input_Reader inp;
	inp.read_file(p_model, "Job.inp");

	///////////////////////////////////////////////////////////////////////////////////
	

	////////////////////////////////////////////////////////////////////////////////
	// Assembling Global Stiffness Matrix and Force Vector

	Matrix Global_K = p_model->assemble_global_stiffness_matrix();
	Matrix F = p_model->get_force_vector();

	std::cout << "Global Stiffness Matrix --------------------------------------------" << std::endl;
	Global_K.print_matrix();
	std::cout << "Force Vector -------------------------------------------------------" << std::endl;
	F.print_matrix();

	///////////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////////
	// Solving Model

	Matrix c = p_model->Solve_model(Global_K,F);
	std::cout << "Solution Vector ----------------------------------------------------" << std::endl;
	c.print_matrix();

	/////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////
	// Exporting results 

	p_model->export_results(c,"Job.out");
	std::cout << "Results Exported sucessfully ---------------------------------------" << std::endl;
	////////////////////////////////////////////////////////////////////////////////

	*/

	Matrix sample{ 2,2 };
	sample.data = {2,5,12,5};

	sample.print_matrix();
	std::cout << "---------------------------------" << std::endl;

	Matrix invert = sample.Invert_Jacobain(sample);
	invert.print_matrix();
	std::cout << "---------------------------------" << std::endl;

	Matrix identity = invert * sample;
	identity.print_matrix();
}