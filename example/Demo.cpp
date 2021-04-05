#include <vector>

#include <mio/mmap.hpp>
#include <AUROC.hpp>

#include "EdgeOnlyCore.hpp"
#include "EdgeNodeCore.hpp"

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

	const bool shouldExportRawScore = false;
	const auto alpha = 1;
	const auto beta = 1;
	const auto gamma = 0.5;
	const auto zeta = 0.7;
	const int shapeCMS[] = {2, 3000};

	// Random seed
	// --------------------------------------------------------------------------------

#ifdef NDEBUG
	const unsigned seed = time(nullptr);
	printf("Seed = %u\t// In case of reproduction\n", seed);
	srand(seed);
#else
	printf("// Debug build is deterministic");
#endif

	// Read dataset
	// --------------------------------------------------------------------------------

	std::error_code err;
	const auto fileMeta = mio::make_mmap_source(pathMeta, err);
	if (err) {
		printf("%s:%d fileMeta: %s\n", __FILE__, __LINE__, err.message().c_str());
		exit(EXIT_FAILURE);
	}
	const auto n = atoi(fileMeta.data());

	const auto fileData = mio::make_mmap_source(pathData, err);
	if (err) {
		printf("%s:%d fileData: %s\n", __FILE__, __LINE__, err.message().c_str());
		exit(EXIT_FAILURE);
	}
	const auto fileLabel = mio::make_mmap_source(pathLabel, err);
	if (err) {
		printf("%s:%d fileLabel: %s\n", __FILE__, __LINE__, err.message().c_str());
		exit(EXIT_FAILURE);
	}
	const auto fileScore = fopen(SOLUTION_DIR"out/Score.tsv", "wb");
	if (!fileScore) {
		printf("%s:%d fileScore: Cannot create file\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

	// Do the magic
	// --------------------------------------------------------------------------------

	int source, destination, timestamp;
	auto itData = fileData.begin() - 1;
	auto itLabel = fileLabel.begin() - 1;
	const auto score = new double[n];
	const auto label = new double[n];
	Isconna::EdgeOnlyCore isc(shapeCMS[0], shapeCMS[1], zeta);
	// Isconna::EdgeNodeCore isc(shapeCMS[0], shapeCMS[1], zeta);
	for (int i = 0; i < n; i++) {
		source = strtol(itData + 1, const_cast<char**>(&itData), 10);
		destination = strtol(itData + 1, const_cast<char**>(&itData), 10);
		timestamp = strtol(itData + 1, const_cast<char**>(&itData), 10);
		double f, w, g;
		isc(source, destination, timestamp, f, w, g);
		score[i] = pow(f, alpha) * pow(w, beta) * pow(g, gamma);
		label[i] = strtol(itLabel + 1, const_cast<char**>(&itLabel), 10);
		if (shouldExportRawScore)
			fprintf(fileScore, "%f\t%f\t%f\t%f\n", score[i], f, w, g);
	}
	fclose(fileScore);
	printf("# Records = %d\t// Process is done\n", n);
	if (shouldExportRawScore)
		printf("// Raw scores are exported to\n// %s\n", SOLUTION_DIR"out/Score.tsv");

	// Evaluate scores
	// --------------------------------------------------------------------------------

	printf("ROC-AUC = %.4f\n", AUROC(label, score, n));
	delete[] score;
	delete[] label;
}
