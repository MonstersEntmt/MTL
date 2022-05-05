#pragma once

#include "ComplexNumber.h"

namespace MTL
{
	struct CubicRoots
	{
	public:
		constexpr CubicRoots() : r0(), r1(), r2() {}
		constexpr CubicRoots(ComplexNumber r0) : r0(r0), r1(), r2() {}
		constexpr CubicRoots(ComplexNumber r0, ComplexNumber r1) : r0(r0), r1(r1), r2() {}
		constexpr CubicRoots(ComplexNumber r0, ComplexNumber r1, ComplexNumber r2) : r0(r0), r1(r1), r2(r2) {}

		ComplexNumber r0, r1, r2;
	};

	// TODO(MarcasRealAccount): Reimplement such that the A variable is equal to
	// { { 4.0, 0.0 }, { -3.7320508075688772935274463415059, 0.0 }, { -0.26794919243112270647255365849413, 0.0 } }
	constexpr CubicRoots SolveCubicRoots(float a, float b, float c, float d)
	{
		float d0 = b * b - 3.0f * a * c;
		float d1 = 2.0f * b * b * b - 9.0f * a * b * c + 27.0f * a * a * d;

		ComplexNumber A = SqrtC(d1 * d1 - 4.0f * d0 * d0 * d0);
		ComplexNumber C = CbrtC((d1 * A) * 0.5f);

		float         factor = -1.0f / (3.0f * a);
		ComplexNumber D      = factor * (b + d0 / C);
		CubicRoots    roots;
		roots.r0 = D + factor * C;
		roots.r1 = D + factor * C * PowC((-1.0f + SqrtC(-3.0f)) * 0.5f, 1);
		roots.r2 = D + factor * C * PowC((-1.0f + SqrtC(-3.0f)) * 0.5f, 2);
		return roots;
	}

	static constexpr CubicRoots A = SolveCubicRoots(1.0f, 0.0f, -15.0f, -4.0f);
} // namespace MTL