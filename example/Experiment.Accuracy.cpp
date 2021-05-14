#include <random>

#include <AUROC.hpp>
#include <mio/mmap.hpp>
#include <sqlite3.h>
#include <tbb/parallel_for_each.h>

#include "EdgeOnlyCore.hpp"
#include "EdgeNodeCore.hpp"

struct Result {
	double expFreq;
	double expWidth;
	double expGap;
	double decay;
	int row;
	int col;
	double auroc;
};

int main(int argc, char* argv[]) {
	// Parameter
	// --------------------------------------------------------------------------------

	// const auto nameDataset = "DARPA";
	// const auto pathMeta = DATASET_DIR"DARPA/processed/Meta.txt";
	// const auto pathData = DATASET_DIR"DARPA/processed/Data.csv";
	// const auto pathLabel = DATASET_DIR"DARPA/processed/Label.csv";

	const auto nameDataset = "CIC-IDS2018";
	const auto pathMeta = DATASET_DIR"CIC-IDS2018/processed/Meta.txt";
	const auto pathData = DATASET_DIR"CIC-IDS2018/processed/Data.csv";
	const auto pathLabel = DATASET_DIR"CIC-IDS2018/processed/Label.csv";

	// const auto nameDataset = "UNSW-NB15";
	// const auto pathMeta = DATASET_DIR"UNSW-NB15/processed/Meta.txt";
	// const auto pathData = DATASET_DIR"UNSW-NB15/processed/Data.csv";
	// const auto pathLabel = DATASET_DIR"UNSW-NB15/processed/Label.csv";

	// const auto nameDataset = "ISCX-IDS2012";
	// const auto pathMeta = DATASET_DIR"ISCX-IDS2012/processed/Meta.txt";
	// const auto pathData = DATASET_DIR"ISCX-IDS2012/processed/Data.csv";
	// const auto pathLabel = DATASET_DIR"ISCX-IDS2012/processed/Label.csv";

	// const auto nameDataset = "CTU-13";
	// const auto pathMeta = DATASET_DIR"CTU-13/processed/Meta.txt";
	// const auto pathData = DATASET_DIR"CTU-13/processed/Data.csv";
	// const auto pathLabel = DATASET_DIR"CTU-13/processed/Label.csv";

	// const auto nameDataset = "CIC-DDoS2019";
	// const auto pathMeta = DATASET_DIR"CIC-DDoS2019/processed/Meta.txt";
	// const auto pathData = DATASET_DIR"CIC-DDoS2019/processed/Data.csv";
	// const auto pathLabel = DATASET_DIR"CIC-DDoS2019/processed/Label.csv";

	const auto expsFreq = {1.0};
	const auto expsWidth = {1.0};
	const auto expsGap = {0.5, 1.0};
	const auto decays = {0.0, 0.3, 0.5, 0.7};
	const std::initializer_list<int[2]> shapesCMS = {{2, 3000}};

	const int numRepeat = 11;

	// Random seed
	// --------------------------------------------------------------------------------

	std::random_device eng;
	std::uniform_int_distribution<> randint;

	// Read meta (total number of records)
	// --------------------------------------------------------------------------------

	std::error_code err;
	const auto fileMeta = mio::make_mmap_source(pathMeta, err);
	if (err) {
		printf("%s:%d fileMeta: %s\n", __FILE__, __LINE__, err.message().c_str());
		exit(EXIT_FAILURE);
	}
	const auto n = atoi(fileMeta.data());

	// Load dataset
	// --------------------------------------------------------------------------------

	const auto src = new int[n];
	const auto dst = new int[n];
	const auto ts = new int[n];
	const auto label = new double[n];
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
	auto it = fileData.begin() - 1;
	for (int i = 0; i < n; i++) {
		src[i] = strtol(it + 1, const_cast<char**>(&it), 10);
		dst[i] = strtol(it + 1, const_cast<char**>(&it), 10);
		ts[i] = strtol(it + 1, const_cast<char**>(&it), 10);
		label[i] = fileLabel[2 * i] == '1';
	}
	printf("# Records = %d\t// Dataset is loaded\n", n);

	// Do the magic
	// --------------------------------------------------------------------------------

	const char* nameAlg;
	std::vector<Result> results;
	tbb::parallel_for_each(expsFreq, [&](double expFreq) {
		tbb::parallel_for_each(expsWidth, [&](double expWidth) {
			tbb::parallel_for_each(expsGap, [&](double expGap) {
				if (expFreq == 0 && expWidth == 0 && expGap == 0) return;
				tbb::parallel_for_each(decays, [&](double decay) {
					tbb::parallel_for_each(shapesCMS, [&](const int shapeCMS[2]) {
						const auto auroc = new double[numRepeat];
						tbb::parallel_for(0, numRepeat, [&](int rep) {
							srand(randint(eng));
							Isconna::EdgeOnlyCore isc(shapeCMS[0], shapeCMS[1], decay);
							// Isconna::EdgeNodeCore isc(shapeCMS[0], shapeCMS[1], decay);
							nameAlg = isc.nameAlg;
							const auto score = new double[n];
							for (int i = 0; i < n; i++)
								score[i] = isc(src[i], dst[i], ts[i], expFreq, expWidth, expGap);
							auroc[rep] = AUROC(label, score, n);
							delete[] score;
						});
						std::sort(auroc, auroc + numRepeat);
						const auto medianAUROC = auroc[numRepeat / 2];
						printf("%g\t%g\t%g\t%g\t%d\t%d\t%.4f\n", expFreq, expWidth, expGap, decay, shapeCMS[0], shapeCMS[1], medianAUROC);
						results.push_back({expFreq, expWidth, expGap, decay, shapeCMS[0], shapeCMS[1], medianAUROC});
						delete[] auroc;
					});
				});
			});
		});
	});

	// Save result
	// --------------------------------------------------------------------------------

	sqlite3* db;
	sqlite3_open_v2(SOLUTION_DIR"out/Result.sqlite", &db, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_WAL, nullptr);
	sqlite3_busy_timeout(db, 1000); // language=sql
	sqlite3_exec(db, "CREATE TABLE IF NOT EXISTS [AUROC.Isconna] (\n\tdataset TEXT,\n\talg TEXT,\n\tf REAL,\n\tw REAL,\n\tg REAL,\n\tdecay REAL,\n\trow INTEGER,\n\tcol INTEGER,\n\troc_auc REAL,\n\tCONSTRAINT [PrimaryKey_AUROC.Isconna]\n\t\tPRIMARY KEY (dataset, alg, f, w, g, decay, row, col)\n);", nullptr, nullptr, nullptr);
	sqlite3_stmt* stmt; // language=sql
	sqlite3_prepare_v2(db, "REPLACE INTO [AUROC.Isconna]\nVALUES (:dataset, :alg, :f, :w, :g, :decay, :row, :col, :roc_auc);", -1, &stmt, nullptr);
	for (const auto& result: results) {
		sqlite3_reset(stmt);
		sqlite3_bind_text(stmt, sqlite3_bind_parameter_index(stmt, ":dataset"), nameDataset, -1, SQLITE_STATIC);
		sqlite3_bind_text(stmt, sqlite3_bind_parameter_index(stmt, ":alg"), nameAlg, -1, SQLITE_STATIC);
		sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, ":f"), result.expFreq);
		sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, ":w"), result.expWidth);
		sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, ":g"), result.expGap);
		sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, ":decay"), result.decay);
		sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, ":row"), result.row);
		sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, ":col"), result.col);
		sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, ":roc_auc"), result.auroc);
		if (sqlite3_step(stmt) != SQLITE_DONE)
			printf("%s:%d %s\n", __FILE__, __LINE__, sqlite3_errmsg(db));
	}
	sqlite3_finalize(stmt); // language=sql
	sqlite3_exec(db, "UPDATE [AUROC.Isconna]\nSET alg='E'\nWHERE alg = 'Isconna-EO';\n\nUPDATE [AUROC.Isconna]\nSET alg='EN'\nWHERE alg = 'Isconna-EN';", nullptr, nullptr, nullptr);
	sqlite3_close(db);

	// Clean up
	// --------------------------------------------------------------------------------

	delete[] src;
	delete[] dst;
	delete[] ts;
	delete[] label;
}
