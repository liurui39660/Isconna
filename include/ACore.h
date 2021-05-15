#pragma once

namespace Isconna {
struct ACore { // A for abstract
	virtual double operator()(int src, int dst, int ts, double alpha, double beta, double gamma) = 0;
	virtual ~ACore() = default; // I know I can move shared methods here, but after all, this class is only for the convenience of dynamic cast.
};
}
