#pragma once

#include <iostream>
#include <vector>
#include <iomanip>

class Matrix
{
public:
	// variables
	int rows{0};
	int cols{0};
	std::vector <double> data;
	
	Matrix() = delete;
	Matrix(int r, int c);


	double& operator()(int i, int j);

	Matrix operator+(Matrix& other) ;
	Matrix operator*(Matrix& other);

	Matrix Transpose();
	void print_matrix()const;
	void PLUDecomposition (Matrix& P , Matrix& L , Matrix& U);


};
