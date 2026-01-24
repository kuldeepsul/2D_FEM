#pragma once
#include "Material.h"

Material::Material(int material_id_param, double E_param, double nu_param)
{
	this->M_id = material_id_param;
	this->E = E_param;
	this->nu = nu_param;

}


