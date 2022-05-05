#pragma once

#include "Math.h"

#include <ostream>
#include <sstream>
#include <type_traits>

namespace MTL
{
	struct ComplexNumber
	{
	public:
		constexpr friend ComplexNumber SinC(ComplexNumber z);
		constexpr friend ComplexNumber CosC(ComplexNumber z);
		constexpr friend ComplexNumber TanC(ComplexNumber z);
		constexpr friend ComplexNumber CotC(ComplexNumber z);
		constexpr friend ComplexNumber SecC(ComplexNumber z);
		constexpr friend ComplexNumber CscC(ComplexNumber z);

		constexpr friend ComplexNumber SinHC(ComplexNumber z);
		constexpr friend ComplexNumber CosHC(ComplexNumber z);
		constexpr friend ComplexNumber TanHC(ComplexNumber z);
		constexpr friend ComplexNumber CotHC(ComplexNumber z);
		constexpr friend ComplexNumber SecHC(ComplexNumber z);
		constexpr friend ComplexNumber CscHC(ComplexNumber z);

		constexpr friend ComplexNumber ArcSinC(ComplexNumber z);
		constexpr friend ComplexNumber ArcCosC(ComplexNumber z);
		constexpr friend ComplexNumber ArcTanC(ComplexNumber z);
		constexpr friend ComplexNumber ArcCotC(ComplexNumber z);
		constexpr friend ComplexNumber ArcSecC(ComplexNumber z);
		constexpr friend ComplexNumber ArcCscC(ComplexNumber z);

		constexpr friend ComplexNumber ExpC(ComplexNumber);
		constexpr friend ComplexNumber LogC(ComplexNumber);
		constexpr friend ComplexNumber PowC(ComplexNumber, ComplexNumber);
		constexpr friend ComplexNumber SqrtC(ComplexNumber);
		constexpr friend ComplexNumber CbrtC(ComplexNumber);

		static constexpr bool IsReal(ComplexNumber complex)
		{
			return complex.m_B == 0.0f;
		}

	public:
		constexpr ComplexNumber() : m_A(0.0f), m_B(0.0f) {}
		constexpr ComplexNumber(float a) : m_A(a), m_B(0.0f) {}
		constexpr ComplexNumber(float a, float b) : m_A(a), m_B(b) {}
		constexpr ComplexNumber(float a, float b, float c) : m_A(a), m_B(b) {}

		constexpr ComplexNumber& setA(float a)
		{
			m_A = a;
			return *this;
		}
		constexpr ComplexNumber& setB(float b)
		{
			m_B = b;
			return *this;
		}

		constexpr float a() const { return m_A; }
		constexpr float b() const { return m_B; }

		constexpr float r() const { return Sqrt(m_A * m_A + m_B * m_B); }
		constexpr float phi() const { return ArcTan2(m_B, m_A); }

		constexpr ComplexNumber& operator+=(ComplexNumber other)
		{
			m_A += other.m_A;
			m_B += other.m_B;
			return *this;
		}
		constexpr ComplexNumber& operator-=(ComplexNumber other)
		{
			m_A -= other.m_A;
			m_B -= other.m_B;
			return *this;
		}
		constexpr ComplexNumber& operator*=(float other)
		{
			m_A *= other;
			m_B *= other;
			return *this;
		}
		constexpr ComplexNumber& operator*=(ComplexNumber other)
		{
			float x = m_A;
			float y = m_B;
			m_A     = x * other.m_A - y * other.m_B;
			m_B     = x * other.m_B + y * other.m_A;
			return *this;
		}
		constexpr ComplexNumber& operator/=(float other)
		{
			float invDen = 1.0f / (other * other);
			m_A *= other * invDen;
			m_B *= other * invDen;
			return *this;
		}
		constexpr ComplexNumber& operator/=(ComplexNumber other)
		{
			float x      = m_A;
			float y      = m_B;
			float invDen = 1.0f / (other.m_A * other.m_A + other.m_B * other.m_B);
			m_A          = (x * other.m_A + y * other.m_B) * invDen;
			m_B          = (y * other.m_A + x * other.m_B) * invDen;
			return *this;
		}
		constexpr friend ComplexNumber operator+(float lhs, ComplexNumber rhs) { return { lhs + rhs.m_A, rhs.m_B }; }
		constexpr friend ComplexNumber operator+(ComplexNumber lhs, float rhs) { return { lhs.m_A + rhs, lhs.m_B }; }
		constexpr friend ComplexNumber operator+(ComplexNumber lhs, ComplexNumber rhs) { return { lhs.m_A + rhs.m_A, lhs.m_B + rhs.m_B }; }
		constexpr friend ComplexNumber operator-(float lhs, ComplexNumber rhs) { return { lhs - rhs.m_A, rhs.m_B }; }
		constexpr friend ComplexNumber operator-(ComplexNumber lhs, float rhs) { return { lhs.m_A - rhs, lhs.m_B }; }
		constexpr friend ComplexNumber operator-(ComplexNumber lhs, ComplexNumber rhs) { return { lhs.m_A - rhs.m_A, lhs.m_B - rhs.m_B }; }
		constexpr friend ComplexNumber operator*(float lhs, ComplexNumber rhs) { return { lhs * rhs.m_A, lhs * rhs.m_B }; }
		constexpr friend ComplexNumber operator*(ComplexNumber lhs, float rhs) { return { lhs.m_A * rhs, lhs.m_B * rhs }; }
		constexpr friend ComplexNumber operator*(ComplexNumber lhs, ComplexNumber rhs) { return { lhs.m_A * rhs.m_A - lhs.m_B * rhs.m_B, lhs.m_A * rhs.m_B + lhs.m_B * rhs.m_A }; }
		constexpr friend ComplexNumber operator/(float lhs, ComplexNumber rhs)
		{
			float invDen = 1.0f / (rhs.m_A * rhs.m_A + rhs.m_B * rhs.m_B);
			return { lhs * rhs.m_A * invDen, lhs * rhs.m_B * invDen };
		}
		constexpr friend ComplexNumber operator/(ComplexNumber lhs, float rhs)
		{
			float invDen = 1.0f / (rhs * rhs);
			return { lhs.m_A * rhs * invDen, lhs.m_B * rhs * invDen };
		}
		constexpr friend ComplexNumber operator/(ComplexNumber lhs, ComplexNumber rhs)
		{
			float invDen = 1.0f / (rhs.m_A * rhs.m_A + rhs.m_B * rhs.m_B);
			return { (lhs.m_A * rhs.m_A + lhs.m_B * rhs.m_B) * invDen, (lhs.m_B * rhs.m_A + lhs.m_A * rhs.m_B) * invDen };
		}
		constexpr friend ComplexNumber operator-(ComplexNumber rhs) { return { -rhs.m_A, -rhs.m_B }; }

		friend std::ostream& operator<<(std::ostream& lhs, ComplexNumber rhs)
		{
			std::ostringstream str;
			str << rhs.m_A << " + " << rhs.m_B << 'i';
			return lhs << str.str();
		}

	private:
		float m_A;
		float m_B;
	};

	static constexpr ComplexNumber Imaginary = { 0.0f, 1.0f };

	constexpr ComplexNumber SinC(ComplexNumber z) { return { Sin(z.m_A) * CosH(z.m_B), Cos(z.m_A) * SinH(z.m_B) }; }
	constexpr ComplexNumber CosC(ComplexNumber z) { return { Cos(z.m_A) * CosH(z.m_B), Sin(z.m_A) * SinH(z.m_B) }; }
	constexpr ComplexNumber TanC(ComplexNumber z) { return SinC(z) / CosC(z); }
	constexpr ComplexNumber CotC(ComplexNumber z) { return 1.0f / TanC(z); }
	constexpr ComplexNumber SecC(ComplexNumber z) { return 1.0f / CosC(z); }
	constexpr ComplexNumber CscC(ComplexNumber z) { return 1.0f / SinC(z); }

	constexpr ComplexNumber SinHC(ComplexNumber z) { return -Imaginary * SinC(Imaginary * z); }
	constexpr ComplexNumber CosHC(ComplexNumber z) { return CosC(Imaginary * z); }
	constexpr ComplexNumber TanHC(ComplexNumber z) { return -Imaginary * TanC(Imaginary * z); }
	constexpr ComplexNumber CotHC(ComplexNumber z) { return Imaginary * CotC(Imaginary * z); }
	constexpr ComplexNumber SecHC(ComplexNumber z) { return SecC(Imaginary * z); }
	constexpr ComplexNumber CscHC(ComplexNumber z) { return Imaginary * CscC(Imaginary * z); }

	constexpr ComplexNumber ArcSinC(ComplexNumber z) { return -Imaginary * LogC(z * Imaginary + SqrtC(1.0f - z * z)); }
	constexpr ComplexNumber ArcCosC(ComplexNumber z) { return -Imaginary * LogC(z + SqrtC(1.0f - z * z)); }
	constexpr ComplexNumber ArcTanC(ComplexNumber z) { return (Imaginary / 2.0f) * LogC((Imaginary + z) / (Imaginary - z)); }
	constexpr ComplexNumber ArcCotC(ComplexNumber z) { return (Imaginary / 2.0f) * LogC((Imaginary - z) / (Imaginary + z)); }
	constexpr ComplexNumber ArcSecC(ComplexNumber z) { return -Imaginary * LogC(1.0f / z + SqrtC(1.0f / (z * z) - 1.0f)); }
	constexpr ComplexNumber ArcCscC(ComplexNumber z) { return -Imaginary * LogC(SqrtC(1.0f - 1.0f / (z * z)) + Imaginary / z); }

	constexpr ComplexNumber ExpC(ComplexNumber z)
	{
		float ex = Exp(z.m_A);
		return { ex * Cos(z.m_B), ex * Sin(z.m_B) };
	}
	constexpr ComplexNumber LogC(ComplexNumber z) { return { Log(z.r()), z.phi() }; }
	constexpr ComplexNumber PowC(ComplexNumber b, ComplexNumber p)
	{
		if (ComplexNumber::IsReal(b))
			return ExpC(p * Log(b.m_A));

		ComplexNumber etr = ExpC(p * Log(b.m_A * b.m_A + b.m_B * b.m_B));
		float         phi = b.phi();
		return etr * (CosC(phi * p) + SinC(phi * p));
	}

	constexpr ComplexNumber SqrtC(ComplexNumber complex)
	{
		if (complex.m_B == 0.0f)
		{
			if (complex.m_A >= 0.0f)
				return Sqrt(complex.m_A);
			else
				return { 0.0f, Sqrt(-complex.m_A) };
		}
		float len = Sqrt(complex.m_A * complex.m_A + complex.m_B * complex.m_B);
		return { Sqrt((complex.m_A + len) * 0.5f), Signum(complex.m_B) * Sqrt((-complex.m_A + len) * 0.5f) };
	}

	constexpr ComplexNumber CbrtC(ComplexNumber complex) { return PowC(complex, 0.33333333f); }

	static constexpr ComplexNumber A = LogC(ComplexNumber { 1, 1 });
} // namespace MTL