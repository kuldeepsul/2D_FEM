#pragma once
#include <iostream>
#include "Matrix.h"


class Material
{
public:
	int M_id{ 0 };
	double E{0};							
	double nu{0};
	Matrix D_stress{3,3};
	Matrix D_strain{3,3};

	Material() = default;
	Material(int material_id_param , double E_param , double nu_param);
};