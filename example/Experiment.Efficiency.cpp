#include <chrono>
#include <vector>

#include <mio/mmap.hpp>
#include <sqlite3.h>

#include "EdgeOnlyCore.hpp"
#include "EdgeNodeCore.hpp"

using namespace std::chrono;

int main(int argc, char* argv[]) {
	// Parameter
	// --------------------------------------------------------------------------------

	// const auto nameDataset = "DARPA";
	// const auto pathMeta = DATASET_DIR"DARPA/processed/Meta.txt";
	// const auto pathData = DATASET_DIR"DARPA/processed/Data.csv";

	const auto nameDataset = "CIC-IDS2018";
	const auto pathMeta = DATASET_DIR"CIC-IDS2018/processed/Meta.txt";
	const auto pathData = DATASET_DIR"CIC-IDS2018/processed/Data.csv";

	// const auto nameDataset = "UNSW-NB15";
	// const auto pathMeta = DATASET_DIR"UNSW-NB15/processed/Meta.txt";
	// const auto pathData = DATASET_DIR"UNSW-NB15/processed/Data.csv";

	// const auto nameDataset = "ISCX-IDS2012";
	// const auto pathMeta = DATASET_DIR"ISCX-IDS2012/processed/Meta.txt";
	// const auto pathData = DATASET_DIR"ISCX-IDS2012/processed/Data.csv";

	// const auto nameDataset = "CTU-13";
	// const auto pathMeta = DATASET_DIR"CTU-13/processed/Meta.txt";
	// const auto pathData = DATASET_DIR"CTU-13/processed/Data.csv";

	// const auto nameDataset = "CIC-DDoS2019";
	// const auto pathMeta = DATASET_DIR"CIC-DDoS2019/processed/Meta.txt";
	// const auto pathData = DATASET_DIR"CIC-DDoS2019/processed/Data.csv";

	const auto nameAlg = "Isconna-EO"; // Isconna-EO or Isconna-EN
	const auto numRepeat = 11;

	// Read best parameter
	// --------------------------------------------------------------------------------

	sqlite3* db;
	sqlite3_open_v2(SOLUTION_DIR"out/Result.sqlite", &db, SQLITE_OPEN_READWRITE | SQLITE_OPEN_WAL, nullptr);
	sqlite3_busy_timeout(db, 1000); // language=sql
	sqlite3_exec(db, "CREATE TABLE IF NOT EXISTS [Time.Isconna] (\n\tdataset TEXT,\n\talg TEXT,\n\tf REAL,\n\tw REAL,\n\tg REAL,\n\tdecay REAL,\n\trow INTEGER,\n\tcol INTEGER,\n\texe_time REAL,\n\tCONSTRAINT [PrimaryKey_Time.Isconna]\n\t\tPRIMARY KEY (dataset, alg, f, w, g, decay, row, col)\n);", nullptr, nullptr, nullptr);
	sqlite3_stmt* stmt; // language=sql
	sqlite3_prepare_v2(db, "SELECT f, w, g, decay, row, col\nFROM [AUROC.Isconna]\nWHERE alg = :alg\n  AND dataset = :dataset\nORDER BY roc_auc DESC\nLIMIT 1;", -1, &stmt, nullptr);
	sqlite3_bind_text(stmt, sqlite3_bind_parameter_index(stmt, ":dataset"), nameDataset, -1, SQLITE_STATIC);
	sqlite3_bind_text(stmt, sqlite3_bind_parameter_index(stmt, ":alg"), nameAlg, -1, SQLITE_STATIC);
	if (sqlite3_step(stmt) != SQLITE_ROW) {
		fprintf(stderr, "%s:%d %s\n", __FILE__, __LINE__, sqlite3_errmsg(db));
		exit(EXIT_FAILURE);
	}
	const auto expFreq = sqlite3_column_double(stmt, 0);
	const auto expWidth = sqlite3_column_double(stmt, 1);
	const auto expGap = sqlite3_column_double(stmt, 2);
	const auto decay = sqlite3_column_double(stmt, 3);
	const int shapeCMS[2] = {sqlite3_column_int(stmt, 4), sqlite3_column_int(stmt, 5)};
	sqlite3_finalize(stmt);

	// Read meta (total number of records)
	// --------------------------------------------------------------------------------

	std::error_code err;
	const auto fileMeta = mio::make_mmap_source(pathMeta, err);
	const auto n = atoi(fileMeta.data());

	// Load dataset
	// --------------------------------------------------------------------------------

	const auto src = new int[n];
	const auto dst = new int[n];
	const auto ts = new int[n];
	const auto fileData = mio::make_mmap_source(pathData, err);
	auto it = fileData.begin() - 1;
	for (int i = 0; i < n; i++) {
		src[i] = strtol(it + 1, const_cast<char**>(&it), 10);
		dst[i] = strtol(it + 1, const_cast<char**>(&it), 10);
		ts[i] = strtol(it + 1, const_cast<char**>(&it), 10);
	}
	printf("# Records = %d\t// Dataset is loaded\n", n);

	// Run
	// --------------------------------------------------------------------------------

	const auto score = new double[n];
	const auto times = new time_t[numRepeat];
	for (int rep = 0; rep < numRepeat; rep++) {
		srand(time(nullptr));
		Isconna::ACore* isc;
		if (!strcmp(nameAlg, Isconna::EdgeOnlyCore::nameAlg)) {
			isc = new Isconna::EdgeOnlyCore(shapeCMS[0], shapeCMS[1], decay);
		} else if (!strcmp(nameAlg, Isconna::EdgeNodeCore::nameAlg)) {
			isc = new Isconna::EdgeNodeCore(shapeCMS[0], shapeCMS[1], decay);
		} // else SelfDestruction();
		const auto timeBegin = high_resolution_clock::now();
		for (int i = 0; i < n; i++)
			score[i] = (*isc)(src[i], dst[i], ts[i], expFreq, expWidth, expGap);
		times[rep] = duration_cast<milliseconds>(high_resolution_clock::now() - timeBegin).count();
		printf("%02d %lld\n", rep, times[rep]);
		delete isc;
	}
	std::sort(times, times + numRepeat);

	// Save results
	// --------------------------------------------------------------------------------

	// language=sql
	sqlite3_prepare_v2(db, "REPLACE INTO [Time.Isconna]\nVALUES (:dataset, :alg, :f, :w, :g, :decay, :row, :col, :exe_time);", -1, &stmt, nullptr);
	sqlite3_bind_text(stmt, sqlite3_bind_parameter_index(stmt, ":dataset"), nameDataset, -1, SQLITE_STATIC);
	sqlite3_bind_text(stmt, sqlite3_bind_parameter_index(stmt, ":alg"), nameAlg, -1, SQLITE_STATIC);
	sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, ":f"), expFreq);
	sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, ":w"), expWidth);
	sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, ":g"), expGap);
	sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, ":decay"), decay);
	sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, ":row"), shapeCMS[0]);
	sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, ":col"), shapeCMS[1]);
	sqlite3_bind_int64(stmt, sqlite3_bind_parameter_index(stmt, ":exe_time"), times[numRepeat / 2]);
	if (sqlite3_step(stmt) != SQLITE_DONE) {
		fprintf(stderr, "%s:%d %s\n", __FILE__, __LINE__, sqlite3_errmsg(db));
		exit(EXIT_FAILURE);
	}
	sqlite3_finalize(stmt);
	sqlite3_close(db);

	// Clean up
	// --------------------------------------------------------------------------------

	delete[] score;
	delete[] times;
}
