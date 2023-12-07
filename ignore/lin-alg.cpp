#include "lin-alg.h"

bool isClose(double scalar,double to, double tolerance)
{
	return std::abs(scalar - to) < tolerance;
}

Vector &Vector::plus(const Vector &other)
{
	if (_value.size() != other._value.size()) throw std::exception();
	for (int i = 0; i < _value.size(); i++)
	{
		_value[i] += other._value[i];
	}
	return *this;
}

Vector &Vector::plus(const Vector &&other)
{
	if (_value.size() != other._value.size()) throw std::exception();
	for (int i = 0; i < _value.size(); i++)
	{
		_value[i] += other._value[i];
	}
	return *this;
}

Vector &Vector::minus(const Vector &other)
{
	if (_value.size() != other._value.size()) throw std::exception();
	for (int i = 0; i < _value.size(); i++)
	{
		_value[i] -= other._value[i];
	}
	return *this;
}

Vector &Vector::minus(const Vector &&other)
{
	if (_value.size() != other._value.size()) throw std::exception();
	for (int i = 0; i < _value.size(); i++)
	{
		_value[i] -= other._value[i];
	}
	return *this;
}

double Vector::norm() const
{
	double result = 0;
	std::for_each(_value.begin(), _value.end(), [&result](auto &elem) { result += elem * elem;});
	return std::sqrt(result);
}

Vector &Vector::operator+(const Vector &other)
{
	if (_value.size() != other._value.size()) throw std::exception();
	for (int i = 0; i < _value.size(); i++)
	{
		_value[i] += other._value[i];
	}
	return *this;
}

Vector &Vector::operator+(const Vector &&other)
{
	if (_value.size() != other._value.size()) throw std::exception();
	for (int i = 0; i < _value.size(); i++)
	{
		_value[i] += other._value[i];
	}
	return *this;
}

Vector &Vector::operator-(const Vector &other)
{
	if (_value.size() != other._value.size()) throw std::exception();
	for (int i = 0; i < _value.size(); i++)
	{
		_value[i] -= other._value[i];
	}
	return *this;
}

Vector &Vector::operator-(const Vector &&other)
{
	if (_value.size() != other._value.size()) throw std::exception();
	for (int i = 0; i < _value.size(); i++)
	{
		_value[i] -= other._value[i];
	}
	return *this;
}

Vector &Vector::operator/(const double scalar)
{
	if (isClose(scalar, 0)) throw std::exception();
	for (double & elem : _value)
	{
		elem /= scalar;
	}
	return *this;
}

Vector &Vector::operator*(const Matrix &other)
{

	return {};
}

Vector &Vector::operator*(const Matrix &&other)
{
	return {};
}

Vector &Vector::multiplication(const Matrix &other)
{
	return {};
}

Vector &Vector::multiplication(const Matrix &&other)
{
	return {};
}

Matrix &Matrix::multiplication(const Matrix &other)
{
	return {};
}

Matrix &Matrix::multiplication(const Matrix &&other)
{
	return <#initializer#>;
}

Matrix &Matrix::multiplication(const Vector &other)
{
	return <#initializer#>;
}

Matrix &Matrix::multiplication(const Vector &&other)
{
	return <#initializer#>;
}

Matrix &Matrix::plus(const Matrix &other)
{
	return <#initializer#>;
}

Matrix &Matrix::plus(const Matrix &&other)
{
	return <#initializer#>;
}

Matrix &Matrix::minus(const Matrix &other)
{
	return <#initializer#>;
}

Matrix &Matrix::minus(const Matrix &&other)
{
	return <#initializer#>;
}

double Matrix::norm() const
{
	return 0;
}

Matrix &Matrix::operator+(const Matrix &other)
{
	return <#initializer#>;
}

Matrix &Matrix::operator+(const Matrix &&other)
{
	return <#initializer#>;
}

Matrix &Matrix::operator-(const Matrix &other)
{
	return <#initializer#>;
}

Matrix &Matrix::operator-(const Matrix &&other)
{
	return <#initializer#>;
}

Matrix &Matrix::operator/(double scalar)
{
	return <#initializer#>;
}

Matrix &Matrix::operator*(const Matrix &other)
{
	return <#initializer#>;
}

Matrix &Matrix::operator*(const Matrix &&other)
{
	return <#initializer#>;
}

Matrix &Matrix::operator*(const Vector &other)
{
	return <#initializer#>;
}

Matrix &Matrix::operator*(const Vector &&other)
{
	return <#initializer#>;
}

Matrix &Matrix::transpose()
{
	return <#initializer#>;
}

Matrix &Matrix::inverse()
{
	return <#initializer#>;
}
