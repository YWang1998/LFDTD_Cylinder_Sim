/*

Header file that contains all the global types and parameters declaration

*/

#pragma once

#define WRITE_PROBE 0

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <map>
#include <cstdio>
#include <iomanip>
#include <chrono>
#include <fstream>
#include <memory>
#define _USE_MATH_DEFINES
#include <math.h>
#include <boost/multi_array.hpp>
#include "algorithm"

// Intel MKL 
#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"

#include "Matrix.h"

// Create a 3D array that is nx x ny x nz
typedef boost::multi_array<int, 4> node_array_type; // nodeNum 4D matrix
typedef boost::multi_array<double, 3> double_array_type; // double 3D matrix

enum Solver
{
    _PARDISO = 1,
};


extern const int nx, ny, nz, nx_TFSF, ny_TFSF, nz_TFSF, qstop, numProbe;
extern const double dt;
extern const int tStep;
extern double s, fc, td, tc;
extern double eps0, mu0, v0;
extern int PML, Nnode, NNZ;
extern int pulseType;
extern double R0, Rs, Rl;
extern double scale;
extern double pi;

extern int maxit;
extern double tol;

extern int jCellIndex[6]; // jCellIndex is assumed to be 1-indexed
extern int jCellIndex_Z[6];
extern int jDirecIndex[3];

extern int mCellIndex[6]; // mCellIndex is assumed to be 1-indexed
extern int mCellIndex_Z[6];
extern int mDirecIndex[3];

extern int probeDirecIndex[3];
extern int jResistorIndex[6]; // jResistorIndex is assumed to be 1-indexed
extern int rDirecIndex[3];
extern int probeCell[];