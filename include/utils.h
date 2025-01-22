#ifndef UTILS_H
#define UTILS_H

#include <cstddef>
#include <string>
#include <vector>
#include <array>
#include <iostream>

/**
 * \brief Retrieves the amount of available (free) physical memory on the current system.
 * \return Size in bytes of the free physical memory.
 */
size_t getAvailableMemory();

/**
 * \brief Reads a file of precomputed binary combinations for a given dimensionality.
 *
 * The file used depends on \p dims (from 2D to 9D).
 * Each line of the file contains a binary string (e.g., "0101") representing a combination.
 *
 * \param dims Integer from 2 to 9, inclusive.
 * \return A vector of binary strings representing combinations.
 */
std::vector<std::string> readCombinations(const int& dims);

/**
 * \brief Checks if a given MBR (array of [min,max] intervals in each dimension)
 *        is considered "valid" with respect to a user-defined hyperplane.
 *
 * \param mbr        The multi-dimensional bounding region, size = dims.
 * \param Comb       A list of precomputed binary strings used to test MBR vertices.
 * \param dims       The number of dimensions (excluding the last "virtual" coordinate).
 * \param queryPlane An array of floats describing the hyperplane (e.g., q1 + q2 + ... < 1).
 *
 * \return True if all vertices of the MBR lie below the hyperplane, false otherwise.
 */
bool MbrIsValid(const std::vector<std::array<float, 2>>& mbr,
                const std::vector<std::string>& Comb,
                int dims,
                const float queryPlane[]);

#endif // UTILS_H
