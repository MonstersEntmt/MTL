#pragma once

#include <bit>

namespace MTL
{
	template <class T>
	struct FloatBits;

	template <>
	struct FloatBits<float>
	{
	public:
		static constexpr std::uint32_t MantissaWidth = 23U;
		static constexpr std::uint32_t ExponentWidth = 8U;
		static constexpr std::uint32_t MantissaMask  = (1U << MantissaWidth) - 1U;
		static constexpr std::uint32_t SignMask      = 1U << (ExponentWidth + MantissaWidth);
		static constexpr std::uint32_t ExponentMask  = ~(SignMask | MantissaMask);
		static constexpr std::uint32_t QuietNaNMask  = 0x0040'0000U;

		static constexpr std::uint32_t ExponentBias = (1U << (ExponentWidth - 1U)) - 1U;
		static constexpr std::uint32_t MaxExponent  = (1U << ExponentWidth) - 1U;

		static constexpr std::uint32_t MinSubNormal = 1U;
		static constexpr std::uint32_t MaxSubNormal = (1U << MantissaWidth) - 1U;
		static constexpr std::uint32_t MinNormal    = 1U << MantissaWidth;
		static constexpr std::uint32_t MaxNormal    = ((MaxExponent - 1U) << MantissaWidth) | MaxSubNormal;

		static constexpr float Zero        = 0.0f;
		static constexpr float NegZero     = -1.0f;
		static constexpr float Infinite    = std::numeric_limits<float>::infinity();
		static constexpr float NegInfinite = -Infinite;

		static constexpr float BuildNaN(std::uint32_t v)
		{
			FloatBits bits = Infinite;
			bits.setMantissa(v);
			return bits;
		}

	public:
		constexpr FloatBits(float value) : m_Value(value) {}
		constexpr FloatBits(std::uint32_t value) : m_Value(std::bit_cast<float>(value)) {}

		constexpr operator float() const { return m_Value; }
		constexpr operator std::uint32_t() const { return std::bit_cast<std::uint32_t>(m_Value); }

		constexpr float         value() const { return m_Value; }
		constexpr std::uint32_t uintValue() const { return std::bit_cast<std::uint32_t>(m_Value); }

		constexpr void setSign(bool sign) { m_Value = std::bit_cast<float>((uintValue() & ~SignMask) | (sign << 31U)); }
		constexpr void setMantissa(std::uint32_t mantissa)
		{
			mantissa &= MantissaMask;
			m_Value = std::bit_cast<float>((uintValue() & ~MantissaMask) | mantissa);
		}
		constexpr void setUnbiasedExponent(std::uint16_t exponent)
		{
			exponent = (exponent << MantissaWidth) & ExponentMask;
			m_Value  = std::bit_cast<float>((uintValue() & ~ExponentMask) | exponent);
		}

		constexpr bool          getSign() const { return (uintValue() & SignMask) >> 31U; }
		constexpr std::uint32_t getMantissa() const { return uintValue() & MantissaMask; }
		constexpr std::uint16_t getUnbiasedExponent() const { return static_cast<std::uint16_t>((uintValue() & ExponentMask) >> MantissaWidth); }
		constexpr int           getExponent() const { return static_cast<int>(getUnbiasedExponent()) - ExponentBias; }

		constexpr bool isZero() const { return getMantissa() == 0U && getUnbiasedExponent() == 0U; }
		constexpr bool isInf() const { return getMantissa() == 0U && getUnbiasedExponent() == MaxExponent; }
		constexpr bool isNaN() const { return getMantissa() != 0U && getUnbiasedExponent() == MaxExponent; }
		constexpr bool isInfOrNaN() const { return getUnbiasedExponent() == MaxExponent; }

	private:
		float m_Value;
	};
} // namespace MTL