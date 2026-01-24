#pragma once
#include "Matrix.h"

Matrix::Matrix(int r, int c)
{
	rows = r;
	cols = c;
	data = std::vector <double>(r * c, 0.0);

}

double& Matrix::operator()(int i, int j)
{
	return this->data[i * cols + j ];
}

Matrix Matrix::operator+(Matrix& other) 
{
	Matrix C(this->rows, this->cols);

	for (int i{ 0 }; i < this->rows; i++)
	{
		for (int j{ 0 }; j < this->cols; j++)
		{
			C(i, j) = (*this)(i, j) + other(i, j);
		}
	}
	return C;
}


Matrix Matrix::Transpose()
{
	Matrix T(this->rows ,this->cols);

	for (int i{ 0 }; i < this->rows; i++)
	{
		for (int j{ 0 }; j < this->cols; j++)
		{
			T(j, i) = (*this)(i, j);
		}
	}
	return T;
}

Matrix Matrix::operator*(Matrix& other)
{
	Matrix C(this->rows, this->cols);
	double sum{ 0 };
	for (int i{ 0 }; i < this->rows; i++)
	{
		for (int j{ 0 }; j < other.cols; j++)
		{
			for (int k{ 0 }; k < other.rows; k++)
			{
				sum = sum + ((*this)(i, k) * other(k, j));
			}
			
			C(i, j) = sum;
			sum = 0;
		}
		
	}
	return C;
}

void Matrix::print_matrix() const 
{
	int count{ 0 };
	std::cout << std::right;
	for (int i{ 0 }; i < size(this->data); i++)
	{
		std::cout << std::setw(8)<< std::setprecision(2) << this->data[i];
		count++;
		if (count == this->cols)
		{
			std::cout << std::endl;
			count = 0;
		}
	}
};

void Matrix::PLUDecomposition (Matrix& P, Matrix& L, Matrix& U)
{
	// Create P , L Matrix as identity

	for (int i{ 0 }; i < this->cols; i++)
	{
		for (int k{ 0 }; k < this->rows; k++)
		{
			if (i == k)
			{
				P(i, k) = 1;
				L(i, k) = 1;
			}
			if (i != k)
			{
				P(i, k) = 0;
				L(i, k) = 0;
			}
		}
	}

	U.data = this->data;

	double pivot{ 0 };

	for (int i{ 0 }; i < U.cols; i++)
	{
		if (U(i, i) != 0)
		{
			pivot = U(i, i);
		}
		else
		{
			for (int k = (i + 1); k < U.rows; k++)
			{
				if (U(k, i) != 0)
				{
					pivot = U(k, i);

					// flipping rows
					for (int j{0}; j < U.cols; j++)
					{
						double temp{0};

						temp = U(i, j);
						U(i, j) = U(k, j);
						U(k, j) = temp;

						temp = P(i, k);
						P(i, j) = P(k, j);
						P(k, j) = temp;

					}
					break;
				}
				else
				{
					if (k == U.rows - 1)
					{
						std::cout << "No pivot found for column : " << i << std::endl;
					}
				}

			}
		}

		// Assuming we have a non zero pivot here

		for (int m = i + 1; m < U.rows; m++)
		{
			double val{ 0 };
			val = U(m, i) / U(i, i);
			L(m, i) = val;

			for (int n = i; n < U.cols; n++)
			{
				U(m, n) = U(m, n) - ((val)*U(i, n));
				double tol = 1e-12;
				if (std::abs(U(m, n)) < tol)
				{
					U(m, n) = 0;
				}
			}
		}
	}
}