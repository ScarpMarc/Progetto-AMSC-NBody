#pragma once
#include <array>
#include <iostream>
#include <cmath>
#include <istream>
#include <ostream>
#include "Constants.h"

template <unsigned int dim>
class Vector
{
public:
	/// <summary>
	/// Default constructor: all components set to 0
	/// </summary>
	/// <typeparam name="dim"></typeparam>
	Vector<dim>() : Vector<dim>( {0} ) {}

	/// <summary>
	/// Standard constructor: all components must be provided
	/// </summary>
	/// <typeparam name="dim"></typeparam>
	Vector<dim>(const std::array<double, dim>& components) : components(components) { }

	/// <summary>
	/// Subscript operator.
	/// Retrieves the vector component corresponding to the argument, starting from 0.
	/// </summary>
	/// <typeparam name="dim">The component to retrieve</typeparam>
	double& operator[](const unsigned int& comp)
	{
		return components[comp];
	}

	/// <summary>
	/// Subscript operator (const version).
	/// Retrieves the const vector component corresponding to the argument, starting from 0.
	/// </summary>
	/// <typeparam name="dim">The component to retrieve</typeparam>
	double operator[](const unsigned int& comp) const
	{
		return components[comp];
	}

	/// <summary>
	/// Vector-to-vector sum
	/// </summary>
	Vector<dim>& operator+=(const Vector<dim>& other)
	{
		for (unsigned int i = 0; i < dim; ++i) components[i] += other[i];
		return *this;
	}

	/// <summary>
	/// Vector-to-vector sum
	/// </summary>
	constexpr friend Vector<dim> operator+(const Vector<dim> lhs, const Vector<dim>& rhs)
	{
		// Apparently declaring the first one with no reference is more optimised
		Vector<dim> out(lhs);
		out += rhs;
		return out;
	}

	/// <summary>
	/// Vector-to-vector subtraction
	/// </summary>
	Vector<dim>& operator-=(const Vector<dim>& other)
	{
		for (unsigned int i = 0; i < dim; ++i) components[i] -= other[i];
		return *this;
	}

	/// <summary>
	/// Vector-to-vector subtraction
	/// </summary>
	constexpr friend Vector<dim> operator-(const Vector<dim> lhs, const Vector<dim>& rhs)
	{
		// Apparently declaring the first one with no reference is more optimised
		Vector<dim> out(lhs);
		out -= rhs;
		return out;
	}

	/// <summary>
	/// Dot (scalar) product
	/// </summary>
	constexpr friend double operator*(Vector<dim> lhs, const Vector<dim>& rhs)
	{
		double sum = 0;
		for (unsigned int i = 0; i < dim; ++i) sum += (lhs[i]) * (rhs[i]);
		return sum;
	}

	/// <summary>
	/// Cross (vector) product
	/// </summary>
	constexpr Vector<dim>& operator^=(const Vector<dim>& other)
	{
		// Not using a cascading switch statement in fear that the compiler may actually perform this evaluation...
		if constexpr (dim == 1)
		{
			components[0] = 0;
		}
		else if constexpr  (dim == 2)
		{
			components[2U] = components[0U] * other[1U] - components[1U] * other[0U];
			components[0U] = 0;
			components[1U] = 0;
		}
		else if constexpr  (dim == 3)
		{
			components[0U] = components[1U] * other[2U] - components[2U] * other[1U];
			components[1U] = components[2U] * other[0U] - components[0U] * other[2U];
			components[2U] = components[0U] * other[1U] - components[1U] * other[0U];
		}
		return *this;
	}

	/// <summary>
	/// Cross (vector) product
	/// </summary>
	const friend Vector<dim>& operator^(const Vector<dim>& lhs, const Vector<dim>& rhs)
	{
		if (dim == 1)
		{
			return Vector<1>();
		}
		else if (dim == 2 || dim == 3)
		{
			Vector<dim> out(lhs);
			out ^= rhs;
			return out;
		}
	}

	/// <summary>
	/// Product by a scalar
	/// </summary>
	constexpr Vector<dim>& operator*=(const double& val)
	{
		for (unsigned int i = 0; i < dim; ++i) components[i] *= val;
		return *this;
	}

	/// <summary>
	/// Product by a scalar
	/// </summary>
	constexpr Vector<dim> operator*(const double& val) const
	{
		Vector<dim> out(components);
		out *= val;
		return out;
	}

	/// <summary>
	/// Division by a scalar
	/// </summary>
	constexpr Vector<dim>& operator/=(const double& val)
	{
		for (unsigned int i = 0; i < dim; ++i) components[i] /= val;
		return *this;
	}

	/// <summary>
	/// Division by a scalar
	/// </summary>
	constexpr Vector<dim> operator/(const double& val) const
	{
		Vector<dim> out(components);
		out /= val;
		return out;
	}

	/// <summary>
	/// Inversion operator
	/// </summary>
	/// <returns>A vector with all components flipped in sign</returns>
	constexpr Vector<dim> operator-() const
	{
		std::array<double, dim> new_components;
		for (unsigned int i = 0; i < dim; ++i)
		{
			new_components[i] = -components[i];
		}

		return Vector<dim>(new_components);
	}

	friend std::ostream& operator<<(std::ostream& str, const Vector<dim>& vec)
	{
		for (unsigned int i = 0; i < dim; ++i)
		{
			double this_dim = vec[i];
			str.write(reinterpret_cast<char*>(&this_dim), sizeof(double));
		}
		return str;
	}

	friend std::istream& operator>>(std::istream& str, Vector<dim>& vec)
	{
		// read obj from stream
		for (unsigned int i = 0; i < dim; ++i)
		{
			str.read(reinterpret_cast<char*>(&(vec[i])), sizeof(double));
		}
		return str;
	}

	/// <summary>
	/// Writes the Vector's contents to a buffer.
	/// We assume the buffer has enough space to accomodate <code>dim</code> values.
	/// </summary>
	/// <param name="out_vec">Array to write to</param>
	/// <param name="start">Where to start in the array</param>
	constexpr void assign_to(double* out_vec, const unsigned long start = 0UL) const
	{
		for (unsigned int i = 0; i < dim; ++i)
		{
			out_vec[i + start] = components[dim];
		}
	}

	/// <summary>
	/// Writes the Vector's contents to a buffer of <code>float</code>s.
	/// We assume the buffer has enough space to accomodate <code>dim</code> values.
	/// </summary>
	/// <param name="out_vec">Array to write to</param>
	/// <param name="start">Where to start in the array</param>
	constexpr void assign_to_float(float* out_vec, const unsigned long start = 0UL) const
	{
		for (unsigned int i = 0; i < dim; ++i)
		{
			out_vec[i + start] = static_cast<float>(components[dim]);
		}
	}

	/// <summary>
	/// Euclidean norm of the vector
	/// </summary>
	double euNorm() const;

private:
	std::array<double, dim> components;
};

#ifndef _WIN32 // VS supports only OpenMP 2.0 and a few OMP3.0 directives. This is unsupported.
#pragma omp declare reduction (VectorSum : Vector<DIM> : omp_out = omp_out + omp_in) initializer (omp_priv=Vector<DIM>())
#endif

template <unsigned int dim>
std::ostream& operator<<(std::ostream& out, const Vector<dim>& vec)
{
	out << "[ ";
	for (unsigned int i = 0; i < dim; ++i) out << vec[i] << " ";
	out << "]";

	return out;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MEMBER FUNCTIONS
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <unsigned int dim>
double Vector<dim>::euNorm() const
{
	double sum = 0;
	for (unsigned int i = 0; i < dim; ++i) sum += components[i] * components[i];
	return std::sqrt(sum);
}