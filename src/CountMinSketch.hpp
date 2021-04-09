#pragma once

#include <cstring>

namespace Isconna {

const int primes3k[] = {3037, 3181, 3307, 3413, 3539, 3643, 3769, 3907};
const int primes7k[] = {6997, 7151, 7297, 7459, 7561, 7687, 7829, 7951};

template<class T>
struct CountMinSketch {
	// Fields
	// --------------------------------------------------------------------------------

	const unsigned int m = 347; // Yes, a magic number, I just pick a random prime
	const int r, c, len;
	int* const param;
	T* const data;

	// Methods
	// --------------------------------------------------------------------------------

	CountMinSketch() = delete;
	CountMinSketch(const CountMinSketch&) = delete;
	CountMinSketch& operator=(const CountMinSketch&) = delete;

	CountMinSketch(int row, int col):
		r(row), c(col), len(r * c),
		param(new int[2 * r]),
		data(new T[len]) {
		for (int i = 0; i < r; i++) {
#ifdef NDEBUG
			param[i] = rand() + 1; // ×0 is not a good idea, see Hash()
			param[r + i] = rand();
#else
			param[i] = primes3k[i];
			param[r + i] = primes7k[i];
#endif
		}
		memset(data, 0, len * sizeof(T));
	}

	CountMinSketch(int numRow, int numColumn, const int* param):
		r(numRow), c(numColumn), len(r * c),
		param(new int[2 * r]),
		data(new T[len]) {
		memcpy(this->param, param, 2 * r * sizeof(int));
		memset(data, 0, len * sizeof(T));
	}

	~CountMinSketch() {
		delete[] param;
		delete[] data;
	}

	void CopyFrom(const CountMinSketch<T>& other) const {
		memcpy(data, other.data, other.len * sizeof(T));
	}

	void ClearAll() const {
		memset(data, 0, len * sizeof(T));
	}

	void MultiplyAll(T by) const {
		for (int i = 0, I = len; i < I; i++) // Vectorization
			data[i] *= by;
	}

	void Hash(int* indexOut, int a, int b = 0) const {
		for (int i = 0; i < r; i++)
			indexOut[i] = i * c + ((a + m * b) * param[i] + param[r + i]) % c;
	}

	T operator()(const int* index) const {
		T least = data[index[0]];
		for (int i = 1; i < r; i++)
			if (least > data[index[i]])
				least = data[index[i]];
		return least;
	}

	T& operator[](int index) {
		return data[index];
	}

	int ArgMin(const int* index) const {
		int indexOut = index[0];
		T least = data[index[0]];
		for (int i = 1; i < r; i++)
			if (least > data[index[i]]) {
				least = data[index[i]];
				indexOut = index[i];
			}
		return indexOut;
	}

	void Add(const int* index, T by = T{1}) const {
		for (int i = 0; i < r; i++)
			data[index[i]] += by;
	}
};
}
