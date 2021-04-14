#pragma once

#include <cmath>
#include <cstring> // memset

namespace Isconna {
struct EdgeOnlyCore {
	static constexpr char nameAlg[] = "Isconna-EO";
	const int row, col; // The same-layout assumption makes this convenient
	const double zeta; // Scale factor, how much to keep
	int tsInternal = 1; // Check if a new timestamp comes
	unsigned* const index; // Due to the same-layout assumption, I only hash once per edge record
	unsigned* const param; // Hashing parameters, just in case it overflows
	bool* bCur, * bAcc; // Busy indicators, use the bool type
	double* fCur, * fAcc, * wCur, * wAcc, * gCur, * gAcc;
	int* wTime, * gTime; // Not timestamps, I explained this in the paper

	EdgeOnlyCore(int row, int col, double zeta = 0):
		row(row), col(col), zeta(zeta),
		index(new unsigned[row]), param(new unsigned[2 * row]),
		bCur(new bool[row * col]), bAcc(new bool[row * col]),
		fCur(new double[row * col]), fAcc(new double[row * col]),
		wCur(new double[row * col]), wAcc(new double[row * col]), wTime(new int[row * col]),
		gCur(new double[row * col]), gAcc(new double[row * col]), gTime(new int[row * col]) {
		for (int i = 0; i < row; i++) {
			param[i] = rand() + 1; // An unfortunate 0 will index all objects to the same cell
			param[i + row] = rand();
		}
		for (int i = 0, I = row * col; i < I; i++) {
			bCur[i] = bAcc[i] = false;
			fCur[i] = fAcc[i] = wCur[i] = wAcc[i] = gCur[i] = gAcc[i] = 0;
			wTime[i] = gTime[i] = 1;
		}
	}

	virtual ~EdgeOnlyCore() {
		delete[] index;
		delete[] param;
		delete[] bCur;
		delete[] bAcc;
		delete[] fCur;
		delete[] fAcc;
		delete[] wCur;
		delete[] wAcc;
		delete[] gCur;
		delete[] gAcc;
		delete[] wTime;
		delete[] gTime;
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

	double operator()(int src, int dst, int ts, double alpha, double beta, double gamma) {
		if (tsInternal < ts) {
			for (int i = 0, I = row * col; i < I; i++) // Vectorization
				fCur[i] *= zeta;
			for (int i = 0, I = row * col; i < I; i++) { // No vectorization
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
			tsInternal = ts;
		}
		for (int i = 0; i < row; i++) {
			index[i] = i * col + ((src + 347 * dst) * param[i] + param[i + row]) % col; // Hashing, 347 is a magic number
			fCur[index[i]]++; // CMS add
			fAcc[index[i]]++; // CMS add
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
		return pow(GTest(Query(fCur), Query(fAcc), ts), alpha) * pow(GTest(wCur[wIndex], wAcc[wIndex], wTime[wIndex]), beta) * pow(GTest(gCur[gIndex], gAcc[gIndex], gTime[gIndex]), gamma);
	}
};
}
