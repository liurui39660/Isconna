#pragma once

namespace Isconna {
struct ACore { // A for abstract
	struct CMSGroup {
		std::valarray<bool> bCur, bAcc; // Busy indicators, use the bool type
		std::valarray<int> wTime, gTime; // Not timestamps, I explained this in the paper
		std::valarray<double> fCur, fAcc, wCur, wAcc, gCur, gAcc;

		explicit CMSGroup(int len): bCur(len), bAcc(len), wTime(1, len), gTime(1, len), fCur(len), fAcc(len), wCur(len), wAcc(len), gCur(len), gAcc(len) {}
	};

	const int row, col; // The same-layout assumption makes this convenient
	const double zeta; // Scale factor, how much to keep
	int tsInternal = 1; // Check if a new timestamp comes
	std::valarray<unsigned> param; // Hashing parameters, just in case it overflows

	ACore(int row, int col, double zeta): row(row), col(col), zeta(zeta), param(2 * row) {
		for (int i = 0; i < row; i++) {
			param[i] = rand() + 1; // An unfortunate 0 will index all objects to the same cell
			param[i + row] = rand();
		}
	}

	static double GTest(double c, double a, double t) {
		return c == 0 || a == 0 || t <= 1 ? 0 : 2 * c * std::abs(std::log(c * (t - 1) / a));
	}

	virtual void Reset(CMSGroup& cms) const {
		cms.fCur *= zeta;
		for (int i = 0, I = row * col; i < I; i++) { // No vectorization
			if (!cms.bCur[i]) { // Doesn't occur in this timestamp
				if (cms.bAcc[i]) { // But occurred in last time timestamp
					cms.gAcc[i] += cms.gCur[i];
					cms.gCur[i] *= zeta;
					cms.gTime[i]++;
				}
				cms.gCur[i]++; // Gap +1
			}
		}
		std::swap(cms.bAcc, cms.bCur); // Now becomes history
		cms.bCur = false;
	}

	virtual void Update(int a, int b, double& fSc, double& wSc, double& gSc, CMSGroup& cms) const {
		double fMinCur = INFINITY, fMinAcc = INFINITY; // Result of CMS Query
		auto wMinTime = INT_MAX, gMinTime = INT_MAX; // INT_MAX
		unsigned wIndex, gIndex; // Result of CMS ArgQuery
		for (int i = 0; i < row; i++) {
			const auto j = i * col + ((a + 347 * b) * param[i] + param[i + row]) % col; // CMS Hash; 347 is a magic number
			fMinCur = std::min(fMinCur, ++cms.fCur[j]); // CMS Add and CMS Query
			fMinAcc = std::min(fMinAcc, ++cms.fAcc[j]); // CMS Add and CMS Query
			if (!cms.bCur[j]) { // Haven't seen this edge in this timestamp
				if (!cms.bAcc[j]) { // This edge didn't occur in last timestamp
					cms.wAcc[j] += cms.wCur[j];
					cms.wCur[j] *= zeta;
					cms.wTime[j]++;
				}
				cms.wCur[j]++; // Width +1
				cms.bCur[j] = true; // Now it's seen
			}
			if (cms.wTime[j] < wMinTime) { // CMS ArgQuery
				wMinTime = cms.wTime[j];
				wIndex = j;
			}
			if (cms.gTime[j] < gMinTime) { // CMS ArgQuery
				gMinTime = cms.gTime[j];
				gIndex = j;
			}
		}
		fSc = GTest(fMinCur, fMinAcc, tsInternal);
		wSc = GTest(cms.wCur[wIndex], cms.wAcc[wIndex], wMinTime);
		gSc = GTest(cms.gCur[gIndex], cms.gAcc[gIndex], gMinTime);
	}

	virtual double operator()(int src, int dst, int ts, double alpha, double beta, double gamma) = 0;
};
}
