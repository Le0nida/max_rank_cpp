#ifndef CSVUTILS_H
#define CSVUTILS_H

#include "geom.h"

/**
 * Reads a CSV file containing records.
 * @param filename   Path to the CSV file.
 * @param numRecords Expected number of records.
 * @param dimensions Number of dimensions in the dataset.
 * @return Vector of Points read from the file.
 */
std::vector<Point> readCSV(const std::string& filename, int numRecords, int dimensions);


/**
 * Reads the query file (one query per line).
 * @param filename  Path to the file.
 * @param numQueries Expected number of queries.
 * @return Vector of query indices.
 */
std::vector<int> readQuery(const std::string& filename, int numQueries);

/**
 * Writes data to a CSV file, supporting both int and double types.
 * @tparam T         Data type (int or double).
 * @param filename   Path to the output file.
 * @param data       Matrix of data (rows: records, columns: attributes).
 * @param headers    Column headers for the output file.
 */
template <typename T>
void writeCSV(const std::string& filename, const std::vector<std::vector<T>>& data, const std::vector<std::string>& headers);

#endif //CSVUTILS_H
