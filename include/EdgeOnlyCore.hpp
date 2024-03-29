#pragma once

#include <valarray>

#include "ACore.hpp"

namespace Isconna {
struct EdgeOnlyCore: ACore {
	static constexpr char nameAlg[] = "Isconna-EO";
	CMSGroup edge;

	EdgeOnlyCore(int row, int col, double zeta): ACore(row, col, zeta), edge(row * col) {}

	double operator()(int src, int dst, int ts, double alpha, double beta, double gamma) override {
		if (tsInternal < ts) {
			Reset(edge);
			tsInternal = ts;
		}
		double fSc, wSc, gSc;
		Update(src, dst, edge, fSc, wSc, gSc);
		return pow(fSc, alpha) * pow(wSc, beta) * pow(gSc, gamma);
	}
};
}
