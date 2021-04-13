#pragma once

#include <cmath>

#include "CountMinSketch.hpp"

namespace Isconna {
struct EdgeOnlyCore {
	static constexpr char nameAlg[] = "Isconna-EO";
	int t = 1;
	const double zeta;
	int* const index;
	CountMinSketch<bool> bCur, bAcc;
	CountMinSketch<double> fCur, fAcc, wCur, wAcc, gCur, gAcc;
	CountMinSketch<int> wTime, gTime;

	EdgeOnlyCore(int row, int col, double zeta = 0):
		zeta(zeta), index(new int[row]),
		bCur(row, col), bAcc(row, col, bCur.param),
		fCur(row, col, bCur.param), fAcc(row, col, bCur.param),
		wCur(row, col, bCur.param), wAcc(row, col, bCur.param), wTime(row, col, bCur.param),
		gCur(row, col, bCur.param), gAcc(row, col, bCur.param), gTime(row, col, bCur.param) {
		for (int i = 0, I = wTime.len; i < I; i++) {
			wTime[i] = gTime[i] = 1;
		}
	}

	virtual ~EdgeOnlyCore() { delete[] index; }

	static double ComputeScore(double c, double a, double t) {
		return c == 0 || a == 0 || t <= 1 ? 0 : 2 * c * std::abs(std::log(c * (t - 1) / a));
	}

	void operator()(int s, int d, int t, double& fSc, double& wSc, double& gSc) {
		if (this->t < t) {
			fCur.MultiplyAll(zeta);
			for (int i = 0; i < bCur.len; i++) {
				if (!bCur[i]) {
					if (bAcc[i]) {
						gAcc[i] += gCur[i];
						gCur[i] *= zeta;
						gTime[i]++;
					}
					gCur[i]++;
				}
			}
			bAcc.CopyFrom(bCur);
			bCur.ClearAll();
			this->t = t;
		}
		fCur.Hash(index, s, d);
		fCur.Add(index);
		fAcc.Add(index);
		for (int i = 0; i < bCur.r; i++) {
			if (!bCur[index[i]]) {
				if (!bAcc[index[i]]) {
					wAcc[index[i]] += wCur[index[i]];
					wCur[index[i]] *= zeta;
					wTime[index[i]]++;
				}
				wCur[index[i]]++;
				bCur[index[i]] = true;
			}
		}
		const auto wIndex = wTime.ArgMin(index);
		const auto gIndex = gTime.ArgMin(index);
		fSc = ComputeScore(fCur(index), fAcc(index), t);
		wSc = ComputeScore(wCur[wIndex], wAcc[wIndex], wTime[wIndex]);
		gSc = ComputeScore(gCur[gIndex], gAcc[gIndex], gTime[gIndex]);
	}
};
}
