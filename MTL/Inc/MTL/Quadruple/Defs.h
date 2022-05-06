#pragma once

#include <cstdint>

#include <compare>

namespace MTL::Quadruple
{
	struct uint128
	{
	public:
		constexpr uint128() : m_L(0ULL), m_H(0ULL) {}
		constexpr uint128(std::uint64_t value) : m_L(value), m_H(0ULL) {}
		constexpr uint128(std::uint64_t l, std::uint64_t h) : m_L(l), m_H(h) {}

		constexpr uint128& operator+=(uint128 rhs)
		{
			if (rhs.m_L >= ((1ULL << 63ULL) - m_L))
				++m_H;
			m_L += rhs.m_L;
			m_H += rhs.m_H;
			return *this;
		}

		constexpr uint128& operator-=(uint128 rhs)
		{
			if (rhs.m_L > m_L)
				--m_H;
			m_L -= rhs.m_L;
			m_H -= rhs.m_H;
			return *this;
		}

		constexpr uint128& operator*=(uint128 rhs);
		constexpr uint128& operator/=(uint128 rhs);

		constexpr friend uint128 operator+(uint128 lhs, uint128 rhs)
		{
			if (rhs.m_L >= ((1ULL << 63ULL) - lhs.m_L))
				++lhs.m_H;
			lhs.m_L += rhs.m_L;
			lhs.m_H += rhs.m_H;
			return lhs;
		}

		constexpr friend uint128 operator-(uint128 lhs, uint128 rhs)
		{
			if (rhs.m_L > lhs.m_L)
				--lhs.m_H;
			lhs.m_L -= rhs.m_L;
			lhs.m_H -= rhs.m_H;
			return lhs;
		}

		constexpr friend uint128 operator*(uint128 lhs, uint128 rhs);
		constexpr friend uint128 operator/(uint128 lhs, uint128 rhs);

		constexpr friend uint128 operator~(uint128 rhs)
		{
			rhs.m_L = ~rhs.m_L;
			rhs.m_H = ~rhs.m_H;
			return rhs;
		}
		constexpr friend uint128 operator^(uint128 lhs, uint128 rhs)
		{
			lhs.m_L ^= rhs.m_L;
			lhs.m_H ^= rhs.m_H;
			return lhs;
		}
		constexpr friend uint128 operator&(uint128 lhs, uint128 rhs)
		{
			lhs.m_L &= rhs.m_L;
			lhs.m_H &= rhs.m_H;
			return lhs;
		}
		constexpr friend uint128 operator|(uint128 lhs, uint128 rhs)
		{
			lhs.m_L |= rhs.m_L;
			lhs.m_H |= rhs.m_H;
			return lhs;
		}
		constexpr friend uint128 operator<<(uint128 lhs, int offset)
		{
			if (offset < 0)
				return lhs >> (-offset);
			else if (offset == 0)
				return lhs;
			else if (offset >= 128)
				return {};

			if (offset >= 64)
			{
				lhs.m_H = lhs.m_L << (offset - 64);
				lhs.m_L = 0U;
			}
			else
			{
				lhs.m_H <<= offset;
				lhs.m_H |= lhs.m_L >> (64U - offset);
				lhs.m_L <<= offset;
			}
			return lhs;
		}
		constexpr friend uint128 operator>>(uint128 lhs, int offset)
		{
			if (offset < 0)
				return lhs << (-offset);
			else if (offset == 0)
				return lhs;
			else if (offset > 128)
				return {};

			if (offset >= 64)
			{
				lhs.m_L = lhs.m_H >> (offset - 64);
				lhs.m_H = 0U;
			}
			else
			{
				lhs.m_L >>= offset;
				lhs.m_L |= lhs.m_H << (64U - offset);
				lhs.m_H >>= offset;
			}
			return lhs;
		}

	public:
		std::uint64_t m_L;
		std::uint64_t m_H;
	};

	struct quadruple
	{
	public:
		static constexpr std::uint64_t TotalBits           = 128U;
		static constexpr std::uint64_t Precision           = 113U;
		static constexpr std::uint64_t ExponentBits        = 15U;
		static constexpr std::uint64_t MantissaBits        = 112U;
		static constexpr std::uint64_t ExponentBias        = (1U << (ExponentBits - 1U)) - 1U;
		static constexpr std::int64_t  NegMaxExponent      = 1U - ExponentBias;
		static constexpr std::uint64_t MaxExponent         = 0x7FFFU - ExponentBias;
		static constexpr std::uint64_t MaxUnbiasedExponent = (1U << ExponentBits) - 1U;

		// Normal representation:
		// s = signbit
		// E = unbiased exponent
		// bias = exponent bias
		// p = precision
		// T = mantissa
		// v = (-1)^s * 2^(E - bias) * (1 + 2^(1 - p) * T)
		// Subnormal representation:
		// emin = 0 - bias
		// v = (-1)^s * 2^(emin) * (0 + 2^(1 - p) * T)

	public:
		constexpr quadruple() : m_MantissaL(0ULL), m_MantissaH(0ULL), m_Exponent(0ULL), m_Sign(0ULL) {}
		constexpr quadruple(std::int64_t integer) {}
		constexpr quadruple(long double floating) {}
		constexpr quadruple(std::uint64_t mantissaL, std::uint64_t mantissaH, std::uint16_t exponent, bool sign) : m_MantissaL(mantissaL), m_MantissaH(mantissaH), m_Exponent(exponent), m_Sign(sign) {}

		constexpr quadruple& operator+=(quadruple rhs);
		constexpr quadruple& operator-=(quadruple rhs);
		constexpr quadruple& operator/=(quadruple rhs);
		constexpr quadruple& operator*=(quadruple rhs);

		constexpr friend quadruple operator+(quadruple lhs, quadruple rhs);
		constexpr friend quadruple operator-(quadruple lhs, quadruple rhs);
		constexpr friend quadruple operator/(quadruple lhs, quadruple rhs);
		constexpr friend quadruple operator*(quadruple lhs, quadruple rhs);
		constexpr friend quadruple operator-(quadruple rhs)
		{
			rhs.m_Sign = !rhs.m_Sign;
			return rhs;
		}

		constexpr friend std::strong_ordering operator<=>(quadruple lhs, quadruple rhs);

		constexpr void setMantissa(uint128 mantissa)
		{
			m_MantissaL = mantissa.m_L;
			m_MantissaH = mantissa.m_H;
		}
		constexpr uint128       getMantissa() const { return { m_MantissaL, m_MantissaH }; }
		constexpr void          setBiasedExponent(std::uint16_t exponent) { m_Exponent = exponent; }
		constexpr std::uint16_t getBiasedExponent() const { return m_Exponent; }
		constexpr void          setExponent(std::int16_t exponent) { m_Exponent = static_cast<std::uint32_t>(exponent) + ExponentBias; }
		constexpr std::int16_t  getExponent() const { return static_cast<std::int32_t>(m_Exponent) - ExponentBias; }
		constexpr void          setSign(bool sign) { m_Sign = sign; }
		constexpr bool          getSign() const { return m_Sign; }

		constexpr bool isInf() const { return m_Exponent == MaxUnbiasedExponent && m_MantissaL == 0U && m_MantissaH == 0U; }
		constexpr bool isNaN() const { return m_Exponent == MaxUnbiasedExponent && (m_MantissaL != 0U || m_MantissaH != 0U); }
		constexpr bool isInfOrNaN() const { return m_Exponent == MaxUnbiasedExponent; }

	private:
		std::uint64_t m_MantissaL : 64;
		std::uint64_t m_MantissaH : 48;
		std::uint64_t m_Exponent : 15;
		std::uint64_t m_Sign : 1;
	};

	static constexpr quadruple Zero        = { 0ULL, 0ULL, 0U, false };
	static constexpr quadruple NegZero     = { 0ULL, 0ULL, 0U, true };
	static constexpr quadruple Infinity    = { 0ULL, 0ULL, 0x7FFF, false };
	static constexpr quadruple NegInfinity = { 0ULL, 0ULL, 0x7FFF, true };

	static constexpr quadruple Max    = { 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFULL, 0x7FFE, false };
	static constexpr quadruple Min    = { 0ULL, 0ULL, 1U, false };
	static constexpr quadruple NegMax = { 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFULL, 0x7FFE, true };
	static constexpr quadruple NegMin = { 0ULL, 0ULL, 1U, true };

	static constexpr quadruple MaxSubNormal    = { 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFULL, 0U, false };
	static constexpr quadruple MinSubNormal    = { 1ULL, 0ULL, 0U, false };
	static constexpr quadruple NegMaxSubNormal = { 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFULL, 0U, true };
	static constexpr quadruple NegMinSubNormal = { 1ULL, 0ULL, 0U, true };

	constexpr quadruple RoundToIntegralTiesToEven(quadruple x);
	constexpr quadruple RoundToIntegralTiesToAway(quadruple x);
	constexpr quadruple RoundToIntegralTowardZero(quadruple x);
	constexpr quadruple RoundToIntegralTowardPositive(quadruple x);
	constexpr quadruple RoundToIntegralTowardNegative(quadruple x);
	constexpr quadruple RoundToIntegralExact(quadruple x);

	constexpr quadruple NextUp(quadruple x);
	constexpr quadruple NextDown(quadruple x);

	constexpr quadruple Remainder(quadruple x, quadruple y);

	constexpr quadruple MinNum(quadruple x, quadruple y);
	constexpr quadruple MaxNum(quadruple x, quadruple y);
	constexpr quadruple MinNumMag(quadruple x, quadruple y);
	constexpr quadruple MaxNumMag(quadruple x, quadruple y);

	constexpr quadruple FusedMultiplyAdd(quadruple x, quadruple y, quadruple z) { return (x * y) + z; }

	constexpr quadruple Exp(quadruple x);
	constexpr quadruple Log(quadruple y);
	constexpr quadruple Pow(quadruple b, quadruple p);
	constexpr quadruple Sqrt(quadruple x);
} // namespace MTL::Quadruple