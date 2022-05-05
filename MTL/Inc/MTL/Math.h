#pragma once

#include "Float.h"
#include "LUTs/Exp.h"
#include "LUTs/Log.h"
#include "LUTs/Sin.h"

#include <cmath>
#include <cstdint>

#include <bit>

namespace MTL
{
	namespace Details
	{
		constexpr std::uint32_t AbsTop12(float x) { return (FloatBits<float>(x).uintValue() >> 20) & 0x7FF; }

		constexpr void SinCosPoly(double x, double x2, const LUTs::SinCos& p, int n, float& sin, float& cos)
		{
			double x3, x4, x5, x6, s, c, c1, c2, s1;

			x4 = x2 * x2;
			x3 = x2 * x;
			c2 = p.c3 + x2 * p.c4;
			s1 = p.s2 + x2 * p.s3;

			c1 = p.c0 + x2 * p.c1;
			x5 = x3 * x2;
			x6 = x4 * x2;

			s = x + x3 * p.s1;
			c = c1 + x4 * p.c2;

			double a = s + x5 * s1;
			double b = c + x6 * c2;
			if (n & 1)
			{
				cos = a;
				sin = b;
			}
			else
			{
				sin = a;
				cos = b;
			}
		}

		constexpr float SinPoly(double x, double x2, const LUTs::SinCos& p, int n)
		{
			double x3, x4, x6, x7, s, c, c1, c2, s1;

			if ((n & 1) == 0)
			{
				x3 = x * x2;
				s1 = p.s2 + x2 * p.s3;

				x7 = x3 * x2;
				s  = x + x3 * p.s1;

				return s + x7 * s1;
			}
			else
			{
				x4 = x2 * x2;
				c2 = p.c3 + x2 * p.c4;
				c1 = p.c0 + x2 * p.c1;

				x6 = x4 * x2;
				c  = c1 + x4 * p.c2;

				return c + x6 * c2;
			}
		}

		constexpr double ReduceFast(double x, const LUTs::SinCos& p, int& n)
		{
			double r;
			r = x * p.hpi_inv;
			n = (static_cast<std::int32_t>(r) + 0x800000) >> 24;
			return x - n * p.hpi;
		}

		constexpr double ReduceLarge(std::uint32_t xi, int& n)
		{
			const std::uint32_t* arr   = &LUTs::InvPiO4[(xi >> 26) & 15];
			int                  shift = (xi >> 23) & 7;
			std::uint64_t        res0, res1, res2;

			xi = (xi & 0xFFFFFF) | 0x800000;
			xi <<= shift;

			res0 = xi * arr[0];
			res1 = static_cast<std::uint64_t>(xi) * arr[4];
			res2 = static_cast<std::uint64_t>(xi) * arr[8];
			res0 = (res2 >> 32) | (res0 << 32);
			res0 += res1;

			n = static_cast<int>((res0 + (1ULL << 61)) >> 62);
			res0 -= static_cast<std::uint64_t>(n) << 62;
			double x = static_cast<std::int64_t>(res0);
			return x * LUTs::Pi63;
		}

		template <class T>
		constexpr void ForceEval(T x)
		{
			[[unused]] volatile T y = x;
		}

		template <class T>
		constexpr T MultiplyAdd(T x, T y, T z)
		{
			return x * y + z;
		}
	} // namespace Details

	static constexpr double EulersNumber  = 2.718281828459045;
	static constexpr float  EulersNumberF = 2.7182818f;
	static constexpr double TwoPi         = 6.283185307179586;
	static constexpr float  TwoPiF        = 6.28318530f;
	static constexpr double Pi            = 3.1415926535897932;
	static constexpr float  PiF           = 3.1415926f;
	static constexpr double HalfPi        = 1.5707963267948966;
	static constexpr float  HalfPiF       = 1.5707963f;
	static constexpr double QuarterPi     = 0.7853981633974483;
	static constexpr float  QuarterPiF    = 0.7853981f;

	template <class T>
	constexpr T MultiplyAdd(T x, T y, T z)
	{
		return x * y + z;
	}

	template <class T>
	constexpr T PolyEval(T x, T a0)
	{
		return a0;
	}

	template <class T, class... Ts>
	constexpr T PolyEval(T x, T a0, Ts... a)
	{
		return MultiplyAdd(x, PolyEval(x, a...), a0);
	}

	constexpr float Exp(float x)
	{
		if constexpr (!std::is_constant_evaluated())
		{
			return std::expf(x);
		}
		else
		{
			FloatBits<float> xbits = x;

			std::uint32_t x_u   = xbits;
			std::uint32_t x_abs = x_u & 0x7FFF'FFFFU;

			// Exceptional values
			// x = -0x1.6D7B18p+5f
			if (x_u == 0xC236'BD8CU) [[unlikely]]
				return 0x1.108A58p-66f - x * 0x1.0p-95f;

			// When |x| >= 89, |x| < 2^-25, or x is nan
			if (x_abs >= 0x42B2'0000U || x_abs <= 0x3280'0000U) [[unlikely]]
			{
				// |x| < 2^-25
				if (xbits.getUnbiasedExponent() <= 101)
					return 1.0f + x;

				// When x < log(2^-150) or nan
				if (x_u >= 0xC2CF'F1B5U)
				{
					// exp(-Inf) = 0
					if (xbits.isInf())
						return 0.0f;
					// exp(nan) = nan
					if (xbits.isNaN())
						return x;
					if constexpr (std::numeric_limits<float>::round_style == std::round_toward_infinity)
						return FloatBits<float>(FloatBits<float>::MinSubNormal);
					// errno = ERANGE;
					return 0.0f;
				}
				// x >= 89 or nan
				if (!xbits.getSign() && x_u >= 0x42B2'0000)
				{
					// x is finite
					if (x_u < 0x7F80'0000U)
					{
						if constexpr (std::numeric_limits<float>::round_style == std::round_toward_neg_infinity || std::numeric_limits<float>::round_style == std::round_toward_zero)
							return FloatBits<float>::MaxNormal;
						// errno = ERANGE;
					}
					return x + FloatBits<float>::Infinite;
				}
			}

			// For -104 < x < 89, to compute exp(x), we perform the following range
			// reduction: find hi, mid, lo such that:
			//   x = hi + mid + lo, in which
			//     hi is an integer,
			//     mid * 2^7 is an integer
			//     -2^(-8) <= lo < 2^-8.
			// In particular,
			//   hi + mid = round(x * 2^7) * 2^(-7).
			// Then,
			//   exp(x) = exp(hi + mid + lo) = exp(hi) * exp(mid) * exp(lo).
			// We store exp(hi) and exp(mid) in the lookup tables LUTs::ExpM1 and LUTs::ExpM2
			// respectively.  exp(lo) is computed using a degree-4 minimax polynomial
			// generated by Sollya.

			// x_hi = (hi + mid) * 2^7 = round(x * 2^7).
			// The default rounding mode for float-to-int conversion in C++ is
			// round toward-zero. To make it round-to-nearest, we add (-1)^sign(x) * 0.5
			// before conversion.
			int x_hi = static_cast<int>(x * 0x1.0p7f + (xbits.getSign() ? -0.5f : 0.5f));
			// Subtract (hi + mid) from x to get lo.
			x -= static_cast<float>(x_hi) * 0x1.0p-7f;
			double xd = static_cast<double>(x);
			x_hi += 104 << 7;
			// hi = x_hi >> 7
			double exp_hi = LUTs::ExpM1[x_hi >> 7];
			// mid * 2^7 = x_hi & 0x0000'007FU;
			double exp_mid = LUTs::ExpM2[x_hi & 0x7F];
			// Degree-4 minimax polynomial generated by Sollya with the following
			// commands:
			//   > display = hexadecimal;
			//   > Q = fpminimax(expm1(x)/x, 3, [|D...|], [-2^-8, 2^-8]);
			//   > Q;
			double exp_lo = PolyEval(xd, 0x1p0, 0x1.FFFFFFFFFF777p-1, 0x1.000000000071Cp-1, 0x1.555566668E5E7p-3, 0x1.55555555EF243p-5);
			double r      = exp_hi * exp_mid * exp_lo;
			if (r > std::numeric_limits<float>::max())
				return FloatBits<float>::Infinite;
			else if (r < -std::numeric_limits<float>::max())
				return FloatBits<float>::NegInfinite;
			else if (r > 0 && r < std::numeric_limits<float>::min())
				return FloatBits<float>::Zero;
			else if (r < 0 && r > -std::numeric_limits<float>::min())
				return FloatBits<float>::NegZero;
			return static_cast<float>(r);
		}
	}
	constexpr float Log(float x)
	{
		if constexpr (!std::is_constant_evaluated())
		{
			return std::logf(x);
		}
		else
		{
			constexpr double c_Log2 = 0x1.62E42FEFA39EFp-1;
			FloatBits<float> xbits  = x;

			switch (xbits.uintValue())
			{
			case 0x4117'8FEBU: // x = 0x1.2F1FD6p+3f
				if (std::numeric_limits<float>::round_style == std::round_to_nearest)
					return 0x1.1FCBCEp+1f;
				break;
			case 0x4C5D'65A5U: // x = 0x1.BACB4Ap+25f
				if (std::numeric_limits<float>::round_style == std::round_to_nearest)
					return 0x1.1E0696p+4f;
				break;
			case 0x65D8'90D3U: // x = 0x1.B121A6p+76f
				if (std::numeric_limits<float>::round_style == std::round_to_nearest)
					return 0x1.A9A3F2p+5f;
				break;
			case 0x6F31'A8ECU: // x = 0x1.6351D8p+95f
				if (std::numeric_limits<float>::round_style == std::round_to_nearest)
					return 0x1.08B512p+6f;
				break;
			case 0x3F80'0001U: // x = 0x1.000002p+0f
				if (std::numeric_limits<float>::round_style == std::round_toward_infinity)
					return 0x1p-23f;
				return 0x1.FFFFFEp-24f;
			case 0x500F'FB03U: // x = 0x1.1FF606p+33f
				if (std::numeric_limits<float>::round_style != std::round_toward_infinity)
					return 0x1.6FDD34p+4f;
				break;
			case 0x7A17'F30AU: // x = 0x1.2FE614p+117f
				if (std::numeric_limits<float>::round_style != std::round_toward_infinity)
					return 0x1.451436p+6f;
				break;
			case 0x5CD6'9E88U: // x = 0x1.AD3D1p+58f
				if (std::numeric_limits<float>::round_style != std::round_toward_infinity)
					return 0x1.45C146p+5f;
				break;
			}

			int m = 0;

			if (xbits.uintValue() < FloatBits<float>::MinNormal || xbits.uintValue() > FloatBits<float>::MaxNormal)
			{
				if (xbits.isZero())
					return FloatBits<float>::NegInfinite;
				if (xbits.getSign() && !xbits.isNaN())
					return FloatBits<float>::BuildNaN(1U << (FloatBits<float>::MantissaWidth - 1U));
				if (xbits.isInfOrNaN())
					return x;
				// Normalize denormal inputs.
				xbits = xbits.value() * 0x1.0p23f;
				m     = -23;
			}

			m += xbits.getExponent();
			// Set bits to 1.m
			xbits.setUnbiasedExponent(0x7FU);
			int f_index = xbits.getMantissa() >> 16U;

			FloatBits<float> f = xbits.uintValue() & ~0x0000'FFFF;

			double d = xbits.value() - f.value();
			d *= LUTs::OneOverF[f_index];

			double extraFactor = MultiplyAdd<double>(m, c_Log2, LUTs::LogF[f_index]);
			double r           = PolyEval(d, extraFactor, 0x1.FFFFFFFFFFFACp-1, -0x1.FFFFFFFEF9CB2p-2, 0x1.5555513BC679Ap-2, -0x1.FFF4805EA441p-3, 0x1.930180DBDE91Ap-3);
			return static_cast<float>(r);
		}
	}

	constexpr void SinCos(float y, float& sin, float& cos)
	{
		FloatBits<float> ybits = y;
		double           x     = y;
		double           s;
		int              n;

		if (Details::AbsTop12(y) < Details::AbsTop12(QuarterPi))
		{
			double x2 = x * x;

			if (Details::AbsTop12(y) < Details::AbsTop12(FloatBits<float>(0x3980'0000U))) // [[unlikely]]
			{
				if (Details::AbsTop12(y) < Details::AbsTop12(FloatBits<float>(0x80'0000U))) // [[unlikely]]
				{
					// Force underflow for tiny y.
					Details::ForceEval<float>(x2);
				}
				sin = y;
				cos = 1.0f;
				return;
			}

			Details::SinCosPoly(x, x2, LUTs::SinCosFTable[0], 0, sin, cos);
		}
		else if (Details::AbsTop12(y) < Details::AbsTop12(120.0f))
		{
			x = Details::ReduceFast(x, LUTs::SinCosFTable[0], n);

			// Setup the signs for sin and cos.
			s = LUTs::SinCosFTable[0].sign[n & 3];

			Details::SinCosPoly(x * s, x * x, LUTs::SinCosFTable[(n & 2) ? 1 : 0], n, sin, cos);
		}
		else if (Details::AbsTop12(y) < Details::AbsTop12(FloatBits<float>::Infinite)) // [[likely]]
		{
			std::uint32_t xi   = ybits.uintValue();
			bool          sign = ybits.getSign();

			x = Details::ReduceLarge(xi, n);

			// Setup signs for sin and cos - include original sign.
			s = LUTs::SinCosFTable[0].sign[(n + sign) & 3];

			Details::SinCosPoly(x * s, x * x, LUTs::SinCosFTable[((n + sign) & 2) ? 1 : 0], n, sin, cos);
		}
		else
		{
			sin = FloatBits<float>::BuildNaN(FloatBits<float>::QuietNaNMask);
			cos = FloatBits<float>::BuildNaN(FloatBits<float>::QuietNaNMask);
		}
	}
	constexpr float Sin(float y)
	{
		if constexpr (!std::is_constant_evaluated())
		{
			return std::sinf(y);
		}
		else
		{
			FloatBits<float> ybits = y;
			double           x     = y;
			double           s;
			int              n;

			if (Details::AbsTop12(y) < Details::AbsTop12(QuarterPi))
			{
				s = x * x;

				if (Details::AbsTop12(y) < Details::AbsTop12(FloatBits<float>(0x3980'0000U))) // [[unlikely]]
				{
					if (Details::AbsTop12(y) < Details::AbsTop12(FloatBits<float>(0x80'0000U))) // [[unlikely]]
					{
						// Force underflow for tiny y.
						Details::ForceEval<float>(s);
					}
					return y;
				}

				return Details::SinPoly(x, s, LUTs::SinCosFTable[0], 0);
			}
			else if (Details::AbsTop12(y) < Details::AbsTop12(120.0f)) // [[likely]]
			{
				x = Details::ReduceFast(x, LUTs::SinCosFTable[0], n);

				// Setup the signs for sin and cos.
				s = LUTs::SinCosFTable[0].sign[n & 3];

				return Details::SinPoly(x * s, x * x, LUTs::SinCosFTable[(n & 2) ? 1 : 0], n);
			}
			else if (Details::AbsTop12(y) < Details::AbsTop12(FloatBits<float>::Infinite))
			{
				std::uint32_t xi   = ybits.uintValue();
				bool          sign = ybits.getSign();

				x = Details::ReduceLarge(xi, n);

				// Setup signs for sin and cos - include original sign
				s = LUTs::SinCosFTable[0].sign[(n + sign) & 3];

				return Details::SinPoly(x * s, x * x, LUTs::SinCosFTable[((n + sign) & 2) ? 1 : 0], n);
			}

			return FloatBits<float>::BuildNaN(FloatBits<float>::QuietNaNMask);
		}
	}
	constexpr float Cos(float y) { return Sin(y + HalfPi); }
	constexpr float Tan(float y)
	{
		float s, c;
		SinCos(y, s, c);
		if (c >= 0.0f && c <= std::numeric_limits<float>::epsilon())
			return FloatBits<float>::Infinite;
		if (c <= 0.0f && c >= -std::numeric_limits<float>::epsilon())
			return FloatBits<float>::NegInfinite;
		return s / c;
	}
	constexpr float Cot(float y)
	{
		float s, c;
		SinCos(y, s, c);
		if (s >= 0.0f && s <= std::numeric_limits<float>::epsilon())
			return FloatBits<float>::Infinite;
		if (s <= 0.0f && s >= -std::numeric_limits<float>::epsilon())
			return FloatBits<float>::NegInfinite;
		return c / s;
	}
	constexpr float Sec(float y)
	{
		float c = Cos(y);
		if (c >= 0.0f && c <= std::numeric_limits<float>::epsilon())
			return FloatBits<float>::Infinite;
		if (c <= 0.0f && c >= -std::numeric_limits<float>::epsilon())
			return FloatBits<float>::NegInfinite;
		return 1.0f / c;
	}
	constexpr float Csc(float y)
	{
		float s = Sin(y);
		if (s >= 0.0f && s <= std::numeric_limits<float>::epsilon())
			return FloatBits<float>::Infinite;
		if (s <= 0.0f && s >= -std::numeric_limits<float>::epsilon())
			return FloatBits<float>::NegInfinite;
		return 1.0f / s;
	}

	constexpr float SinH(float x) { return (1.0f - Exp(-2.0f * x)) / (2.0f * Exp(-x)); }
	constexpr float CosH(float x) { return (1.0f + Exp(-2.0f * x)) / (2.0f * Exp(-x)); }
	constexpr float TanH(float x)
	{
		float e2x = Exp(2.0f * x);
		return (e2x - 1.0f) / (e2x + 1.0f);
	}
	constexpr float CotH(float x)
	{
		float e2x = Exp(2.0f * x);
		return (e2x + 1.0f) / (e2x - 1.0f);
	}
	constexpr float SecH(float x) { return (2.0f * Exp(x)) / (Exp(2.0f * x) + 1); }
	constexpr float CscH(float x) { return (2.0f * Exp(x)) / (Exp(2.0f * x) - 1); }

	constexpr float ArcSin(float x)
	{
		constexpr float piby2_tail  = 7.5497894159e-08F; /* 0x33a22168 */
		constexpr float hpiby2_head = 7.8539812565e-01F; /* 0x3f490fda */
		constexpr float piby2       = 1.5707963705e+00F; /* 0x3fc90fdb */

		FloatBits<float> xbits  = x;
		std::uint32_t    aux    = xbits.uintValue() & 0x7FFF'FFFFU;
		std::uint32_t    xs     = xbits.uintValue() ^ aux;
		float            spiby2 = FloatBits<float>(xs | FloatBits<float>(piby2).uintValue());
		int              xexp   = xbits.getExponent();
		float            y      = FloatBits<float>(aux);

		// abs(x) >= 0.5
		int transform = xexp >= -1;

		float y2 = y * y;
		float rt = 0.5f * (1.0f - y);
		float r  = transform ? rt : y2;

		// Use a rational approximation for [0.0, 0.5]
		float a = Details::MultiplyAdd(r, Details::MultiplyAdd(r, Details::MultiplyAdd(r, -0.00396137437848476485201154797087F, -0.0133819288943925804214011424456F), -0.0565298683201845211985026327361F), 0.184161606965100694821398249421F);

		float b = Details::MultiplyAdd(r, -0.836411276854206731913362287293F, 1.10496961524520294485512696706F);
		float u = r * (a / b);

		float s  = Sqrt(r);
		float s1 = FloatBits<float>(FloatBits<float>(s).uintValue() & 0xFFFF'0000U);
		float c  = Details::MultiplyAdd(-s1, s1, r) / (s + s1);
		float p  = Details::MultiplyAdd(2.0f * s, u, -Details::MultiplyAdd(c, -2.0f, piby2_tail));
		float q  = Details::MultiplyAdd(s1, -2.0f, hpiby2_head);
		float vt = hpiby2_head - (p - q);
		float v  = Details::MultiplyAdd(y, u, y);
		v        = transform ? vt : v;

		float ret = FloatBits<float>(xs | FloatBits<float>(v).uintValue());
		ret       = aux > 0x3F80'0000U ? FloatBits<float>(0x7FC0'0000U).value() : ret;
		ret       = aux == 0x3F80'0000U ? spiby2 : ret;
		ret       = xexp < -14 ? x : ret;

		return ret;
	}
	constexpr float ArcCos(float x)
	{
		constexpr float piby2      = 1.5707963705e+00F;
		constexpr float pi         = 3.1415926535897933e+00F;
		constexpr float piby2_head = 1.5707963267948965580e+00F;
		constexpr float piby2_tail = 6.12323399573676603587e-17F;

		FloatBits<float> xbits = x;
		std::uint32_t    aux   = xbits.uintValue() & 0x7FFF'FFFFU;
		int              xneg  = xbits.uintValue() != aux;
		int              xexp  = xbits.getExponent();
		float            y     = FloatBits<float>(aux);

		// transform if |x| >= 0.5
		int transform = xexp >= -1;

		float y2 = y * y;
		float yt = 0.5f * (1.0f - y);
		float r  = transform ? yt : y2;

		// Use a rational approximation for [0.0, 0.5]
		float a = Details::MultiplyAdd(r, Details::MultiplyAdd(r, Details::MultiplyAdd(r, -0.00396137437848476485201154797087F, -0.0133819288943925804214011424456F), -0.0565298683201845211985026327361F), 0.184161606965100694821398249421F);

		float b = Details::MultiplyAdd(r, -0.836411276854206731913362287293F, 1.10496961524520294485512696706F);
		float u = r * (a / b);

		float s     = Sqrt(r);
		y           = s;
		float s1    = FloatBits<float>(FloatBits<float>(s).uintValue() & 0xFFFF'0000U);
		float c     = Details::MultiplyAdd(s1, -s1, r) / (s + s1);
		float rettn = Details::MultiplyAdd(s + Details::MultiplyAdd(y, u, -piby2_tail), -2.0f, pi);
		float rettp = 2.0F * (s1 + Details::MultiplyAdd(y, u, c));
		float rett  = xneg ? rettn : rettp;
		float ret   = piby2_head - (x - Details::MultiplyAdd(x, -u, piby2_tail));

		ret = transform ? rett : ret;
		ret = aux > 0x3F80'0000U ? FloatBits<float>(0x7FC0'0000U).value() : ret;
		ret = xbits.uintValue() == 0x3F80'0000U ? 0.0f : ret;
		ret = xbits.uintValue() == 0xBF80'0000U ? pi : ret;
		ret = xexp < -26 ? piby2 : ret;
		return ret;
	}
	constexpr float ArcTan(float x)
	{
		constexpr float piby2 = 1.5707963267948966f; // 0x3ff921fb54442d18

		FloatBits<float> xbits = x;
		std::uint32_t    aux   = xbits.uintValue() & 0x7FFF'FFFFU;
		std::uint32_t    sx    = xbits.uintValue() ^ aux;

		float spiby2 = FloatBits<float>(sx | FloatBits<float>(piby2).uintValue());

		float v = FloatBits<float>(aux);

		// Return for NaN
		float ret = x;

		// 2^26 <= |x| <= Inf => atan(x) is close to piby2
		ret = aux <= 0x7F80'0000U ? spiby2 : ret;

		// Reduce arguments 2^-19 <= |x| < 2^26

		// 39/16 <= x < 2^26
		x       = -1.0f / v;
		float c = 1.57079632679489655800f; // atan(infinity)

		// 19/16 <= x < 39/16
		int   l  = aux < 0x401C'0000U;
		float xx = (v - 1.5f) / Details::MultiplyAdd(v, 1.5f, 1.0f);
		x        = l ? xx : x;
		c        = l ? 9.82793723247329054082e-01f : c; // atan(1.5)

		// 11/16 <= x < 19/16
		l  = aux < 0x3F98'0000U;
		xx = (v - 1.0f) / (1.0f + v);
		x  = l ? xx : x;
		c  = l ? 7.85398163397448278999e-01f : c; // atan(1)

		// 7/16 <= x < 11/16
		l  = aux < 0x3F30'0000U;
		xx = Details::MultiplyAdd(v, 2.0f, -1.0f) / (2.0f + v);
		x  = l ? xx : x;
		c  = l ? 4.63647609000806093515e-01f : c; // atan(0.5)

		// 2^-19 <= x < 7/16
		l = aux < 0x3EE0'0000U;
		x = l ? v : x;
		c = l ? 0.0f : c;

		// Core approximation: Remez(2,2) on [-7/16,7/16]

		float s = x * x;
		float a = Details::MultiplyAdd(s, Details::MultiplyAdd(s, 0.470677934286149214138357545549e-2f, 0.192324546402108583211697690500f), 0.296528598819239217902158651186f);

		float b = Details::MultiplyAdd(s, Details::MultiplyAdd(s, 0.299309699959659728404442796915f, 0.111072499995399550138837673349e1f), 0.889585796862432286486651434570f);

		float q = x * s * (a / b);

		float z  = c - (q - x);
		float zs = FloatBits<float>(sx | FloatBits<float>(z).uintValue());

		ret = aux < 0x4C80'0000U ? zs : ret;

		// |x| < 2^-19
		ret = aux < 0x3600'0000U ? x : ret;
		return ret;
	}
	constexpr float ArcCot(float x) { return HalfPi - ArcTan(x); }
	constexpr float ArcSec(float x)
	{
		if (x == 0.0f)
			return FloatBits<float>::Infinite;
		return ArcCos(1.0f / x);
	}
	constexpr float ArcCsc(float x) { return HalfPi - ArcSec(x); }
	constexpr float ArcTan2(float y, float x)
	{
		if (x > 0.0f)
			return ArcTan(y / x);
		else if (y >= 0.0f && x < 0.0f)
			return ArcTan(y / x) + Pi;
		else if (y < 0.0f && x < 0.0f)
			return ArcTan(y / x) - Pi;
		else if (y > 0.0f && x == 0.0f)
			return HalfPi;
		else if (y < 0.0f && x == 0.0f)
			return -HalfPi;
		else
			return FloatBits<float>::BuildNaN(FloatBits<float>::QuietNaNMask);
	}

	constexpr float ArcSinH(float x) { return Log(x + Sqrt(x * x + 1.0f)); }
	constexpr float ArcCosH(float x) { return Log(x + Sqrt(x * x - 1.0f)); }
	constexpr float ArcTanH(float x) { return 0.5f * Log((1.0f + x) / (1.0f - x)); }
	constexpr float ArcCotH(float x) { return 0.5f * Log((1.0f - x) / (1.0f + x)); }
	constexpr float ArcSecH(float x) { return Log((1.0f + Sqrt(1.0f - x * x)) / x); }
	constexpr float ArcCscH(float x) { return Log(1.0f / x + Sqrt(1.0f / (x * x) + 1.0f)); }

	constexpr float Pow(float b, float p) { return Exp(Log(b) * p); }
	constexpr float Sqrt(float x) { return Exp(Log(x) * 0.5f); }
	constexpr float Cbrt(float x) { return Exp(Log(x) * 0.33333333f); }
	constexpr float Qbrt(float x) { return Exp(Log(x) * 0.25f); }

	constexpr float Signum(float x) { return x < 0.0f ? -1.0f : (x > 0.0f ? 1.0f : 0.0f); }

	static constexpr float A = ArcSec(1.0f);
} // namespace MTL