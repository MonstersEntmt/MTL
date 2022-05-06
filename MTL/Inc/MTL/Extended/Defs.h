#pragma once

#include "../Utils/Core.h"

namespace MTL::Extended
{
#if BUILD_IS_TOOLSET_MSVC
	struct extended
	{
	public:

	public:
		std::uint64_t m_L, m_H;
	};
#else
	using extended = long double;
#endif
} // namespace MTL::Extended