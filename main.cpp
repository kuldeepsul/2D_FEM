#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include "Matrix.h"
#include "Model.h"
#include "Material.h"

int main()
{
	model model_01(0);
	Material mat_def;
	mat_def.M_id = 1;
	mat_def.E = 2e08;
	model_01.create_node(1,0,0,0);
	model_01.create_node(2, 5, 5, 0);
	model_01.create_node(3, 10, 0, 0);
	//model_01.create_node(4, 15, 0, 0);
	model_01.create_frame_element(1, 1, 2, &mat_def, 0.1, 0.0005);
	model_01.create_frame_element(2, 2, 3, &mat_def, 0.1, 0.0005);
	model_01.create_frame_element(3, 3, 1, &mat_def, 0.1, 0.0005);

	Matrix matrix_01 = model_01.assemble_global_stiffness_matrix();
	matrix_01.print_matrix();
}