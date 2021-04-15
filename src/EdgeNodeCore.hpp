#pragma once

#include <valarray>

namespace Isconna {
struct EdgeNodeCore {
	static constexpr char nameAlg[] = "Isconna-EN";
	const int row, col;
	const double zeta;
	int tsInternal = 1;
	unsigned* const index;
	unsigned* const param;
	std::valarray<bool> ebCur, ebAcc, sbCur, sbAcc, dbCur, dbAcc;
	std::valarray<double> efCur, efAcc, ewCur, ewAcc, egCur, egAcc;
	std::valarray<double> sfCur, sfAcc, swCur, swAcc, sgCur, sgAcc;
	std::valarray<double> dfCur, dfAcc, dwCur, dwAcc, dgCur, dgAcc;
	std::valarray<int> ewTime, egTime, swTime, sgTime, dwTime, dgTime;

	EdgeNodeCore(int row, int col, double zeta = 0):
		row(row), col(col), zeta(zeta), index(new unsigned[row]), param(new unsigned[2 * row]),
		ebCur(row * col), ebAcc(row * col), sbCur(row * col), sbAcc(row * col), dbCur(row * col), dbAcc(row * col),
		efCur(row * col), efAcc(row * col), ewCur(row * col), ewAcc(row * col), egCur(row * col), egAcc(row * col),
		sfCur(row * col), sfAcc(row * col), swCur(row * col), swAcc(row * col), sgCur(row * col), sgAcc(row * col),
		dfCur(row * col), dfAcc(row * col), dwCur(row * col), dwAcc(row * col), dgCur(row * col), dgAcc(row * col),
		ewTime(1, row * col), egTime(1, row * col), swTime(1, row * col), sgTime(1, row * col), dwTime(1, row * col), dgTime(1, row * col) {
		for (int i = 0; i < row; i++) {
			param[i] = rand() + 1;
			param[i + row] = rand();
		}
	}

	virtual ~EdgeNodeCore() {
		delete[] index;
		delete[] param;
	}

	static double GTest(double c, double a, double t) {
		return c == 0 || a == 0 || t <= 1 ? 0 : 2 * c * std::abs(std::log(c * (t - 1) / a));
	}

	template<class T>
	T Query(const std::valarray<T>& data) const {
		T least = data[index[0]];
		for (int i = 1; i < row; i++)
			least = std::min(least, data[index[i]]);
		return least;
	}

	template<class T>
	unsigned ArgQuery(const std::valarray<T>& data) const {
		unsigned arg = index[0];
		T least = data[arg];
		for (int i = 1; i < row; i++)
			if (least > data[index[i]]) {
				arg = index[i];
				least = data[arg];
			}
		return arg;
	}

	void ResetComponent(
		std::valarray<double>& fCur,
		std::valarray<bool>& bCur, std::valarray<bool>& bAcc,
		std::valarray<double>& gCur, std::valarray<double>& gAcc, std::valarray<int>& gTime
	) const {
		fCur *= zeta;
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
		bCur = false;
	}

	void UpdateComponent(
		int a, int b, double& fSc, double& wSc, double& gSc,
		std::valarray<double>& fCur, std::valarray<double>& fAcc,
		std::valarray<bool>& bCur, const std::valarray<bool>& bAcc,
		std::valarray<double>& wCur, std::valarray<double>& wAcc, std::valarray<int>& wTime,
		std::valarray<double>& gCur, std::valarray<double>& gAcc, std::valarray<int>& gTime
	) const {
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
		return pow(std::max({ewSc, swSc, dwSc}), beta)
			* pow(std::max({efSc, sfSc, dfSc}), alpha)
			* pow(std::max({egSc, sgSc, dgSc}), gamma);
	}
};
}
