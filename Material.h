#pragma once
#include <iostream>


class Material
{
public:
	int M_id{ 0 };
	double E{0};							
	double nu{0};

	Material() = default;
	Material(int material_id_param , double E_param , double nu_param);
};