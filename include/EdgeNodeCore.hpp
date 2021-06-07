#pragma once

#include <valarray>

#include "ACore.hpp"

namespace Isconna {
struct EdgeNodeCore: ACore {
	static constexpr char nameAlg[] = "Isconna-EN";
	CMSGroup edge, source, destination;

	EdgeNodeCore(int row, int col, double zeta): ACore(row, col, zeta), edge(row * col), source(row * col), destination(row * col) {}

	double operator()(int src, int dst, int ts, double alpha, double beta, double gamma) override {
		if (tsInternal < ts) {
			Reset(edge);
			Reset(source);
			Reset(destination);
			tsInternal = ts;
		}
		double efSc, ewSc, egSc, sfSc, swSc, sgSc, dfSc, dwSc, dgSc;
		Update(src, dst, edge, efSc, ewSc, egSc);
		Update(src, 000, source, sfSc, swSc, sgSc);
		Update(dst, 000, destination, dfSc, dwSc, dgSc);
		return pow(std::max({efSc, sfSc, dfSc}), alpha)
			* pow(std::max({ewSc, swSc, dwSc}), beta)
			* pow(std::max({egSc, sgSc, dgSc}), gamma);
	}
};
}
