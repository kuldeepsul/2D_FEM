#pragma once
#include "Material.h"

Material::Material(int material_id_param, double E_param, double nu_param)
{
	this->M_id = material_id_param;
	this->E = E_param;
	this->nu = nu_param;

	double scaler_D_stress = E / (1 - (nu * nu));
	double scaler_D_strain = E / (1 - nu) * (1 - 2 * nu);

	this->D_stress.data = {		1	,	nu	,	0	,
								nu	,	1	,	0	,
								0	,	0	,(1-nu)/2.0 };




	this->D_strain.data = { (1-nu) ,  nu  ,  0 
							, nu   ,(1-nu),  0
							, 0    ,  0   ,  (1-(2*nu))/2.0};

	this->D_strain.scaler_multiple(scaler_D_strain);
	this->D_stress.scaler_multiple(scaler_D_stress);

}


