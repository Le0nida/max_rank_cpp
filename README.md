# Best Rank

This repository contains a C++ implementation of the “Maximum Rank Query” algorithms described by Mouratidis et al. (*Maximum Rank Query*). However, **unlike the original approach** (which focuses on maximization), this version is adapted for **minimization**—hence the name **Best Rank**.

---

## Overview

The **goal** is to compute, for a given set of multi-dimensional records (points) and a set of queries referencing specific records, the “best rank” that each queried record can achieve under **all possible weight vectors** in a minimization scenario. Specifically, the algorithm determines how many other points are “better” or “worse” than the query record across various dimensions by considering the entire space of weight vectors.

---

## Requirements

- **C++17** (or newer) compatible compiler (e.g., GCC, Clang, MSVC).
- [Eigen](http://eigen.tuxfamily.org/) for linear algebra.  
  The project expects Eigen to be installed at:  
  `C:/Program Files (x86)/C-Libraries/eigen-3.4.0`  
  (Adjust the path in the CMake file if necessary.)
- [HiGHS](https://www.highs.dev/) for linear programming routines.  
  The project expects HiGHS to be installed at:  
  `C:/Program Files (x86)/C-Libraries/HiGHS`  
  (Its libraries are located in `C:/Program Files (x86)/C-Libraries/HiGHS/build/bin`.)  
  Adjust the paths in the CMake file if needed.
- A build system such as **CMake** (recommended) to compile the project easily.

---

## Building the Project

1. **Clone or download** this repository.
2. Verify that the paths to Eigen and HiGHS in the `CMakeLists.txt` file match your system's installation. In this project, they are set as follows:
  - **Eigen**: `C:/Program Files (x86)/C-Libraries/eigen-3.4.0`
  - **HiGHS**:
    - Includes: `C:/Program Files (x86)/C-Libraries/HiGHS` (with additional include directories for `src`, `src/lp_data`, `src/util`, and `build`)
    - Library: `C:/Program Files (x86)/C-Libraries/HiGHS/build/bin/libhighs.a`
3. **Configure** the project using CMake.

---

# Usage

## Mandatory Parameters

- **datafile**  
  Path to the CSV file containing your dataset.

- **numRecords**  
  Number of records (rows) to read from the dataset file.

- **dimensions**  
  Number of dimensions in each record (e.g., 2 for 2D, 3 for 3D, etc.).

- **numQueries**  
  Number of queries to read from the query file.

- **queryfile**  
  Path to a file listing query indices (one per line). Each index references a record in the dataset.

- **outdir**  
  Path to an output directory where result files will be saved.

## Optional 7th Argument

- **Config file path** containing lines of the form `key=value`  
  **or**  
  CLI flags in the form `--key=value` (multiple flags allowed).

---

## Example Config File

### my_config.txt

```text
limitHamWeight=999
maxLevelQTree=8
maxCapacityQNode=20
maxNoBinStringToCheck=999999
halfspacesLengthLimit=21
```

---

## Parameters Summary

- **limitHamWeight** (integer, default=999)  
  Restricts combinatorial searches generating fewer combinations of binary strings.

- **maxLevelQTree** (integer, default=99)  
  Maximum allowed depth of the quad-tree structure.

- **maxCapacityQNode** (integer, default=10)  
  Maximum number of halfspaces in a leaf node before splitting further.

- **maxNoBinStringToCheck** (integer, default=999999)  
  Upper bound on how many binary strings to consider in enumerations.

- **halfspacesLengthLimit** (integer, default=21)  
  Restricts combinatorial searches limiting the number of halfspaces to consider in enumerations.

You can pass these either through the config file or via CLI flags. Defaults apply if none are specified.

---

## Input Files Format

### Dataset CSV

A CSV file with exactly **numRecords** rows (plus 1 for headers), each having **index** and **dimensions** columns. For example, a 3D dataset:

```csv
id,x,y,z
1,0.1,0.2,0.3
...
```

### Query File

A text file with **numQueries** lines, each containing one integer index (1-based).

We provide example datasets and query files in the `examples/` directory.

---

## Example Command Lines

**With a Config File (no CLI flags):**

```text
C:\Users\username\Documents\dataset.csv
10000
3
15
C:\Users\username\Documents\queryfile.txt
C:\Users\username\Desktop
C:\Users\username\Documents\config.txt
```
