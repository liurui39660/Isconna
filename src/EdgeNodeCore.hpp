#pragma once

#include <cmath>
#include <algorithm>

namespace Isconna {
struct EdgeNodeCore {
	static constexpr char nameAlg[] = "Isconna-EN";
	const int row, col;
	const double zeta;
	int tsInternal = 1;
	unsigned* const index;
	unsigned* const param;
	bool* ebCur, * ebAcc, * sbCur, * sbAcc, * dbCur, * dbAcc;
	double* efCur, * efAcc, * ewCur, * ewAcc, * egCur, * egAcc;
	double* sfCur, * sfAcc, * swCur, * swAcc, * sgCur, * sgAcc;
	double* dfCur, * dfAcc, * dwCur, * dwAcc, * dgCur, * dgAcc;
	int* ewTime, * egTime, * swTime, * sgTime, * dwTime, * dgTime;

	EdgeNodeCore(int row, int col, double zeta = 0):
		row(row), col(col), zeta(zeta),
		index(new unsigned[row]), param(new unsigned[2 * row]),
		ebCur(new bool[row * col]), ebAcc(new bool[row * col]),
		efCur(new double[row * col]), efAcc(new double[row * col]),
		ewCur(new double[row * col]), ewAcc(new double[row * col]), ewTime(new int[row * col]),
		egCur(new double[row * col]), egAcc(new double[row * col]), egTime(new int[row * col]),
		sbCur(new bool[row * col]), sbAcc(new bool[row * col]),
		sfCur(new double[row * col]), sfAcc(new double[row * col]),
		swCur(new double[row * col]), swAcc(new double[row * col]), swTime(new int[row * col]),
		sgCur(new double[row * col]), sgAcc(new double[row * col]), sgTime(new int[row * col]),
		dbCur(new bool[row * col]), dbAcc(new bool[row * col]),
		dfCur(new double[row * col]), dfAcc(new double[row * col]),
		dwCur(new double[row * col]), dwAcc(new double[row * col]), dwTime(new int[row * col]),
		dgCur(new double[row * col]), dgAcc(new double[row * col]), dgTime(new int[row * col]) {
		for (int i = 0; i < row; i++) {
			param[i] = rand() + 1;
			param[i + row] = rand();
		}
		for (int i = 0, I = row * col; i < I; i++) {
			ebCur[i] = ebAcc[i] = false;
			efCur[i] = efAcc[i] = ewCur[i] = ewAcc[i] = egCur[i] = egAcc[i] = 0;
			ewTime[i] = egTime[i] = 1;
			sbCur[i] = sbAcc[i] = false;
			sfCur[i] = sfAcc[i] = swCur[i] = swAcc[i] = sgCur[i] = sgAcc[i] = 0;
			swTime[i] = sgTime[i] = 1;
			dbCur[i] = dbAcc[i] = false;
			dfCur[i] = dfAcc[i] = dwCur[i] = dwAcc[i] = dgCur[i] = dgAcc[i] = 0;
			dwTime[i] = dgTime[i] = 1;
		}
	}

	virtual ~EdgeNodeCore() {
		delete[] index;
		delete[] param;
		delete[] ebCur;
		delete[] ebAcc;
		delete[] efCur;
		delete[] efAcc;
		delete[] ewCur;
		delete[] ewAcc;
		delete[] ewTime;
		delete[] egCur;
		delete[] egAcc;
		delete[] egTime;
		delete[] sbCur;
		delete[] sbAcc;
		delete[] sfCur;
		delete[] sfAcc;
		delete[] swCur;
		delete[] swAcc;
		delete[] swTime;
		delete[] sgCur;
		delete[] sgAcc;
		delete[] sgTime;
		delete[] dbCur;
		delete[] dbAcc;
		delete[] dfCur;
		delete[] dfAcc;
		delete[] dwCur;
		delete[] dwAcc;
		delete[] dwTime;
		delete[] dgCur;
		delete[] dgAcc;
		delete[] dgTime;
	}

	static double GTest(double c, double a, double t) {
		return c == 0 || a == 0 || t <= 1 ? 0 : 2 * c * std::abs(std::log(c * (t - 1) / a));
	}

	template<class T>
	T Query(const T* data) const {
		T least = data[index[0]];
		for (int i = 1; i < row; i++)
			if (least > data[index[i]])
				least = data[index[i]];
		return least;
	}

	template<class T>
	unsigned ArgQuery(const T* data) const {
		unsigned arg = index[0];
		T least = data[arg];
		for (int i = 1; i < row; i++)
			if (least > data[index[i]]) {
				arg = index[i];
				least = data[arg];
			}
		return arg;
	}

	void ResetComponent(double* fCur, bool*& bCur, bool*& bAcc, double* gCur, double* gAcc, int* gTime) const {
		for (int i = 0, I = row * col; i < I; i++)
			fCur[i] *= zeta;
		for (int i = 0, I = row * col; i < I; i++) {
			if (!bCur[i]) {
				if (bAcc[i]) {
					gAcc[i] += gCur[i];
					gCur[i] *= zeta;
					gTime[i]++;
				}
				gCur[i]++;
			}
		}
		std::swap(bAcc, bCur);
		memset(bCur, 0, row * col * sizeof(bool));
	}

	void UpdateComponent(int a, int b, double& fSc, double& wSc, double& gSc, double* fCur, double* fAcc, bool* bCur, const bool* bAcc, double* wCur, double* wAcc, int* wTime, double* gCur, double* gAcc, int* gTime) const {
		for (int i = 0; i < row; i++) {
			index[i] = i * col + ((a + 347 * b) * param[i] + param[i + row]) % col;
			fCur[index[i]]++;
			fAcc[index[i]]++;
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
		const auto wIndex = ArgQuery(wTime);
		const auto gIndex = ArgQuery(gTime);
		fSc = GTest(Query(fCur), Query(fAcc), tsInternal);
		wSc = GTest(wCur[wIndex], wAcc[wIndex], wTime[wIndex]);
		gSc = GTest(gCur[gIndex], gAcc[gIndex], gTime[gIndex]);
	}

	double operator()(int src, int dst, int ts, double alpha, double beta, double gamma) {
		if (tsInternal < ts) {
			ResetComponent(efCur, ebCur, ebAcc, egCur, egAcc, egTime);
			ResetComponent(sfCur, sbCur, sbAcc, sgCur, sgAcc, sgTime);
			ResetComponent(dfCur, dbCur, dbAcc, dgCur, dgAcc, dgTime);
			tsInternal = ts;
		}
		double efSc, ewSc, egSc, sfSc, swSc, sgSc, dfSc, dwSc, dgSc;
		UpdateComponent(src, dst, efSc, ewSc, egSc, efCur, efAcc, ebCur, ebAcc, ewCur, ewAcc, ewTime, egCur, egAcc, egTime);
		UpdateComponent(src, 000, sfSc, swSc, sgSc, sfCur, sfAcc, sbCur, sbAcc, swCur, swAcc, swTime, sgCur, sgAcc, sgTime);
		UpdateComponent(dst, 000, dfSc, dwSc, dgSc, dfCur, dfAcc, dbCur, dbAcc, dwCur, dwAcc, dwTime, dgCur, dgAcc, dgTime);
		return pow(std::max({efSc, sfSc, dfSc}), alpha) * pow(std::max({ewSc, swSc, dwSc}), beta) * pow(std::max({egSc, sgSc, dgSc}), gamma);
	}
};
}
