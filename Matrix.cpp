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

Matrix Matrix::scaler_multiple(double scaler)
{
	for (int i{ 0 }; i < this->rows; i++)
	{
		for (int j{ 0 }; j < this->cols; j++)
		{
			(*this)(i, j) = scaler * (*this)(i, j);
		}
	}

	return (*this);
}

void Matrix::print_matrix() const 
{
	std::cout << std::right;
	for (int i{ 0 }; i < this->rows; i++)
	{
		for (int j{ 0 }; j < this->cols; j++)
		{
			std::cout << std::setprecision(4) << std::setw(12) << this->data[i*cols + j];
		}
		std::cout << std::endl;
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

void Matrix::remove_row(int row_no)
{
	for (int i{ (this->rows)-1  }; i >= 0; i--)
	{
		for (int j{ (this->cols) -1 }; j >= 0; j--)
		{
			if (j == row_no || i == row_no)
			{
				this->data.erase(this->data.begin() + (i*cols + j));
			}
		}
	}
	this->rows = (this->rows) - 1;
	this->cols = (this->cols) - 1;
}

Matrix Matrix::Solve_linear_system(Matrix& K, Matrix& f,bool up)
{
	double sum{ 0 };
	double val{ 0 };
	double res{ 0 };
	Matrix x(f.rows,1);
	
	if (up)
	{
		for (int i{ 0 }; i < K.rows; i++)
		{
			for (int j{ 0 }; j < K.cols; j++)
			{
				if (j == i)
				{
					val = K(i, j);
				}
				else
				{
					sum = sum + ((K(i, j) )* (x(j, 0)));
				}
			}
			res = (f(i,0) - sum) / val;
			if (std::abs(res) < 1e-12)
			{
				x(i, 0) = 0;
	
			}
			else
			{
				x(i, 0) = res;
				
			}
			sum = 0;
			val = 0;
			
		}
		return x;
	}
	else
	{
		for (int i{K.rows - 1}; i >= 0; i--)
		{
			for (int j{ K.rows - 1}; j >= 0; j--)
			{
				if (j == i)
				{
					val = K(i, j);
				}
				else
				{
					sum = sum + ((K(i, j)) * (x(j, 0)));
				}
			}
			res = (f(i, 0) - sum) / val;
			if (std::abs(res) < 1e-12)
			{
				x(i, 0) = 0;
				
			}
			else
			{
				x(i, 0) = res;
				
			}
			sum = 0;
			val = 0;
			
		}
		return x;
	}
}

Matrix Matrix::Invert_Jacobain(Matrix& Jacobian)
{
	if (Jacobian.rows != 2 || Jacobian.cols != 2)
	{
		std::cout << "Error : Jacobian is not a 2 x 2 matrix " << std::endl;
	}

	double det_j = (Jacobian(0, 0) * Jacobian(1, 1)) - (Jacobian(0, 1) * Jacobian(1, 0));

	if (det_j == 0)
	{
		std::cout << "Jacobian is not invertable." << std::endl;
	}

	Matrix adj_Jacobian{ 2,2 };
	adj_Jacobian.data = { Jacobian(1, 1) ,-Jacobian(0, 1) ,-Jacobian(1, 0), Jacobian(0, 0) };

	return adj_Jacobian.scaler_multiple(1/det_j);

}


