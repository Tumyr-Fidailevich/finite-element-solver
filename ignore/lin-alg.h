#ifndef FINITE_ELEMENT_SOLVER_LIN_ALG_H
#define FINITE_ELEMENT_SOLVER_LIN_ALG_H

#include "pch.h"
#include "exceptions.h"


class Vector;

class Matrix;

bool isClose(double scalar, double to, double tolerance=1e-16);

class Vector
{
	std::vector<double> _value;
public:
	template <typename Iterator>
	Vector(Iterator begin, Iterator end);

	// TODO Определить операторы
	// TODO Определить оператор обращения через квадратные скобки

	Vector &multiplication(const Matrix &other);

	Vector &multiplication(const Matrix &&other);

	Vector &plus(const Vector &other);

	Vector &plus(const Vector &&other);

	Vector &minus(const Vector &other);

	Vector &minus(const Vector &&other);

	[[nodiscard]] double norm() const;


	Vector& operator+(const Vector &other);

	Vector& operator+(const Vector &&other);

	Vector& operator-(const Vector &other);

	Vector& operator-(const Vector &&other);

	Vector& operator/(double scalar);

	Vector& operator*(const Matrix &other);

	Vector& operator*(const Matrix &&other);

	friend Matrix;
};

template<typename Iterator>
Vector::Vector(Iterator begin, Iterator end)
{

}


class Matrix
{
	std::vector<std::vector<double>> _value;

public:

	// TODO Определить метод взятия обратной матрицы
	// TODO Определить операцию транспонирования
	// TODO Определить операторы

	template <typename Iterator>
	Matrix(Iterator begin, Iterator end);

	// TODO Определить операторы
	// TODO Определить оператор обращения через квадратные скобки

	Matrix &multiplication(const Matrix &other);

	Matrix &multiplication(const Matrix &&other);

	Matrix &multiplication(const Vector &other);

	Matrix &multiplication(const Vector &&other);

	Matrix &plus(const Matrix &other);

	Matrix &plus(const Matrix &&other);

	Matrix &minus(const Matrix &other);

	Matrix &minus(const Matrix &&other);

	Matrix &transpose();

	Matrix &inverse();

	[[nodiscard]] double norm() const;


	Matrix& operator+(const Matrix &other);

	Matrix& operator+(const Matrix &&other);

	Matrix& operator-(const Matrix &other);

	Matrix& operator-(const Matrix &&other);

	Matrix& operator/(double scalar);

	Matrix& operator*(const Matrix &other);

	Matrix& operator*(const Matrix &&other);

	Matrix& operator*(const Vector &other);

	Matrix& operator*(const Vector &&other);

	friend Vector;
};

template<typename Iterator>
Matrix::Matrix(Iterator begin, Iterator end)
{

}

#endif