#pragma once
#include <array>

template <size_t dim>
class Vector
{
public:
	/// <summary>
	/// Default constructor: all components set to 0
	/// </summary>
	/// <typeparam name="dim"></typeparam>
	Vector<dim>() : Vector<dim>(double[dim]) {}

	/// <summary>
	/// Standard constructor: all components must be provided
	/// </summary>
	/// <typeparam name="dim"></typeparam>
	Vector<dim>(const std::array<double, dim>& components) : components(components) { }

	/// <summary>
	/// Retrieves the vector component corresponding to the argument, starting from 0.
	/// </summary>
	/// <typeparam name="dim">The component to retrieve</typeparam>
	const double& operator[](const size_t& comp) const
	{
		return components[comp];
	}

	Vector<dim>& operator+=(const Vector<dim>& other)
	{
		for (unsigned int i = 0; i < dim; ++i) components[i] += other[i];
		return this;
	}


	constexpr friend Vector<dim> operator+(const Vector<dim> lhs, const Vector<dim>& rhs)
	{
		// Apparently declaring the first one with no reference is more optimised
		Vector<dim> out(lhs);
		out += rhs;
		return out;
	}

	/// <summary>
	/// Dot (scalar) product
	/// </summary>
	constexpr friend double operator*(Vector<dim> lhs, const Vector<dim>& rhs)
	{
		double sum = 0;
		for (unsigned int i = 0; i < dim; ++i) sum += lhs[i] * rhs[i];
		return sum;
	}

	constexpr Vector<1>& operator^=(const Vector<1>&)
	{
		components[0] = 0;
		return this;
	}

	constexpr friend double operator^(const Vector<1>&, const Vector<1>&)
	{
		return 0;
	}

	constexpr Vector<2>& operator^=(const Vector<2>& other)
	{
		components[2] = this[0] * other[1] - this[1] * other[0];
		components[0] = 0;
		components[1] = 0;
		return this;
	}

	constexpr friend Vector<2> operator^(const Vector<2>& lhs, const Vector<2>& rhs)
	{
		Vector<dim> out(lhs);
		out ^= rhs;
		return out;
	}

	/// <summary>
	/// Cross (vector) product that stores its result into this object.
	/// </summary>
	constexpr Vector<3>& operator^=(const Vector<3>& other)
	{
		components[0] = this[1] * other[2] - this[2] * other[1];
		components[1] = this[2] * other[0] - this[0] * other[2];
		components[2] = this[0] * other[1] - this[1] * other[0];
		return this;
	}

	/// <summary>
	/// Cross (vector) product
	/// </summary>
	constexpr friend Vector<3>& operator^(const Vector<3>& lhs, const Vector<3>& rhs)
	{
		Vector<dim> out(lhs);
		out ^= rhs;
		return out;
	}

	constexpr Vector<3>& operator/(const double& val) const
	{
		std::array<double, dim> out_components = components;
		for (unsigned int i = 0; i < dim; ++i) out_components[i] /= val;
		return Vector<3>(out_components);
	}

	/// <summary>
	/// Euclidean norm of the vector
	/// </summary>
	double eu_norm() const;

private:
	std::array<double, dim> components;
};