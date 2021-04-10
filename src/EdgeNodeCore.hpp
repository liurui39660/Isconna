#pragma once

#include <cmath>
#include <algorithm>

#include "CountMinSketch.hpp"

namespace Isconna {
struct EdgeNodeCore {
	using CMSb = CountMinSketch<bool>;
	using CMSf = CountMinSketch<double>;
	using CMSi = CountMinSketch<int>;
	static constexpr char nameAlg[] = "Isconna-EN";
	int tInternal = 1;
	const double zeta;
	int* const index;
	CMSb ebCur, ebAcc, sbCur, sbAcc, dbCur, dbAcc;
	CMSf efCur, efAcc, ewCur, ewAcc, egCur, egAcc;
	CMSf sfCur, sfAcc, swCur, swAcc, sgCur, sgAcc;
	CMSf dfCur, dfAcc, dwCur, dwAcc, dgCur, dgAcc;
	CMSi ewTime, egTime, swTime, sgTime, dwTime, dgTime;

	EdgeNodeCore(int row, int col, double zeta = 0):
		zeta(zeta), index(new int[row]),
		ebCur(row, col), ebAcc(row, col, ebCur.param),
		sbCur(row, col), sbAcc(row, col, sbCur.param),
		dbCur(row, col), dbAcc(row, col, dbCur.param),
		efCur(row, col, ebCur.param), efAcc(row, col, ebCur.param),
		sfCur(row, col, sbCur.param), sfAcc(row, col, sbCur.param),
		dfCur(row, col, dbCur.param), dfAcc(row, col, dbCur.param),
		ewCur(row, col, ebCur.param), ewAcc(row, col, ebCur.param), ewTime(row, col, ebCur.param),
		egCur(row, col, ebCur.param), egAcc(row, col, ebCur.param), egTime(row, col, ebCur.param),
		swCur(row, col, sbCur.param), swAcc(row, col, sbCur.param), swTime(row, col, sbCur.param),
		sgCur(row, col, sbCur.param), sgAcc(row, col, sbCur.param), sgTime(row, col, sbCur.param),
		dwCur(row, col, dbCur.param), dwAcc(row, col, dbCur.param), dwTime(row, col, dbCur.param),
		dgCur(row, col, dbCur.param), dgAcc(row, col, dbCur.param), dgTime(row, col, dbCur.param) {
		for (int i = 0, I = ewTime.len; i < I; i++) ewTime[i] = egTime[i] = 1;
		for (int i = 0, I = swTime.len; i < I; i++) swTime[i] = sgTime[i] = 1;
		for (int i = 0, I = dwTime.len; i < I; i++) dwTime[i] = dgTime[i] = 1;
	}

	virtual ~EdgeNodeCore() { delete[] index; }

	static double ComputeScore(double c, double a, double t) {
		return c == 0 || a == 0 || t <= 1 ? 0 : 2 * c * std::abs(std::log(c * (t - 1) / a));
	}

	void ResetComponent(
		CMSf& fCur,
		CMSb& bCur, CMSb& bAcc,
		CMSf& gCur, CMSf& gAcc, CMSi& gTime
	) const {
		for (int i = 0; i < fCur.len; i++) {
			fCur[i] *= zeta;
		}
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
	}

	void UpdateComponent(
		int a, int b,
		double& fSc, double& wSc, double& gSc,
		CMSf& fCur, CMSf& fAcc,
		CMSb& bCur, CMSb& bAcc,
		CMSf& wCur, CMSf& wAcc, CMSi& wTime,
		CMSf& gCur, CMSf& gAcc, CMSi& gTime
	) const {
		fCur.Hash(index, a, b);
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
		fSc = ComputeScore(fCur(index), fAcc(index), tInternal);
		wSc = ComputeScore(wCur[wIndex], wAcc[wIndex], wTime[wIndex]);
		gSc = ComputeScore(gCur[gIndex], gAcc[gIndex], gTime[gIndex]);
	}

	void operator()(int s, int d, int t, double& fSc, double& wSc, double& gSc) {
		if (tInternal < t) {
			ResetComponent(efCur, ebCur, ebAcc, egCur, egAcc, egTime);
			ResetComponent(sfCur, sbCur, sbAcc, sgCur, sgAcc, sgTime);
			ResetComponent(dfCur, dbCur, dbAcc, dgCur, dgAcc, dgTime);
			tInternal = t;
		}
		double efSc, ewSc, egSc, sfSc, swSc, sgSc, dfSc, dwSc, dgSc;
		UpdateComponent(s, d, efSc, ewSc, egSc, efCur, efAcc, ebCur, ebAcc, ewCur, ewAcc, ewTime, egCur, egAcc, egTime);
		UpdateComponent(s, 0, sfSc, swSc, sgSc, sfCur, sfAcc, sbCur, sbAcc, swCur, swAcc, swTime, sgCur, sgAcc, sgTime);
		UpdateComponent(d, 0, dfSc, dwSc, dgSc, dfCur, dfAcc, dbCur, dbAcc, dwCur, dwAcc, dwTime, dgCur, dgAcc, dgTime);
		fSc = std::max({efSc, sfSc, dfSc});
		wSc = std::max({ewSc, swSc, dwSc});
		gSc = std::max({egSc, sgSc, dgSc});
	}
};
}
