#include <vector>
#include <chrono>

#include <AUROC.hpp>
#include <mio/mmap.hpp>

#include "EdgeOnlyCore.hpp"
#include "EdgeNodeCore.hpp"

using namespace std::chrono;

int main(int argc, char* argv[]) {
	// Parameter
	// --------------------------------------------------------------------------------

	// const auto pathMeta = DATASET_DIR"DARPA/processed/Meta.txt";
	// const auto pathData = DATASET_DIR"DARPA/processed/Data.csv";
	// const auto pathLabel = DATASET_DIR"DARPA/processed/Label.csv";

	const auto pathMeta = DATASET_DIR"CIC-IDS2018/processed/Meta.txt";
	const auto pathData = DATASET_DIR"CIC-IDS2018/processed/Data.csv";
	const auto pathLabel = DATASET_DIR"CIC-IDS2018/processed/Label.csv";

	// const auto pathMeta = DATASET_DIR"UNSW-NB15/processed/Meta.txt";
	// const auto pathData = DATASET_DIR"UNSW-NB15/processed/Data.csv";
	// const auto pathLabel = DATASET_DIR"UNSW-NB15/processed/Label.csv";

	// const auto pathMeta = DATASET_DIR"ISCX-IDS2012/processed/Meta.txt";
	// const auto pathData = DATASET_DIR"ISCX-IDS2012/processed/Data.csv";
	// const auto pathLabel = DATASET_DIR"ISCX-IDS2012/processed/Label.csv";

	// const auto pathMeta = DATASET_DIR"CTU-13/processed/Meta.txt";
	// const auto pathData = DATASET_DIR"CTU-13/processed/Data.csv";
	// const auto pathLabel = DATASET_DIR"CTU-13/processed/Label.csv";

	// const auto pathMeta = DATASET_DIR"CIC-DDoS2019/processed/Meta.txt";
	// const auto pathData = DATASET_DIR"CIC-DDoS2019/processed/Data.csv";
	// const auto pathLabel = DATASET_DIR"CIC-DDoS2019/processed/Label.csv";

	const auto alpha = 1;
	const auto beta = 1;
	const auto gamma = 0.5;
	const auto zeta = 0.7;
	const int shapeCMS[] = {2, 3000};

	// Random seed
	// --------------------------------------------------------------------------------

	const unsigned seed = time(nullptr);
	srand(seed);
	printf("Seed = %u\t// In case of reproduction\n", seed);

	// Read dataset
	// --------------------------------------------------------------------------------

	std::error_code err;
	const auto fileMeta = mio::make_mmap_source(pathMeta, err);
	const auto fileData = mio::make_mmap_source(pathData, err);
	const auto fileLabel = mio::make_mmap_source(pathLabel, err);
	const auto fileScore = fopen(SOLUTION_DIR"out/Score.txt", "wb");
	const auto n = atoi(fileMeta.data());

	// Do the magic
	// --------------------------------------------------------------------------------

	auto iteratorData = fileData.begin() - 1;
	auto iteratorLabel = fileLabel.begin() - 1;
	const auto score = new double[n];
	const auto label = new double[n];
	Isconna::EdgeOnlyCore isc(shapeCMS[0], shapeCMS[1], zeta);
	// Isconna::EdgeNodeCore isc(shapeCMS[0], shapeCMS[1], zeta);
	printf("# Records = %d\t// Algorithm is started\n", n);
	const auto timeBegin = steady_clock::now();
	for (int i = 0; i < n; i++) {
		const int src = strtol(iteratorData + 1, const_cast<char**>(&iteratorData), 10);
		const int dst = strtol(iteratorData + 1, const_cast<char**>(&iteratorData), 10);
		const int ts = strtol(iteratorData + 1, const_cast<char**>(&iteratorData), 10);
		score[i] = isc(src, dst, ts, alpha, beta, gamma);
		label[i] = static_cast<int>(strtol(iteratorLabel + 1, const_cast<char**>(&iteratorLabel), 10));
		// fprintf(fileScore, "%f\n", score[i]);
	}
	printf("Time = %lldms\t// Process is done\n", duration_cast<milliseconds>((steady_clock::now() - timeBegin)).count());
	if(ftell(fileScore)) // If anything is written
		printf("// Raw scores are exported to\n// %s\n", SOLUTION_DIR"out/Score.txt");
	fclose(fileScore);

	// Evaluate scores
	// --------------------------------------------------------------------------------

	printf("ROC-AUC = %.4f\n", AUROC(label, score, n));
	delete[] score;
	delete[] label;
}
