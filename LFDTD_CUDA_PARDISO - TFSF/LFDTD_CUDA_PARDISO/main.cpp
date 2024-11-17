// includes, system

#include "global.h"
#include "LFDTD_Coe.h"
#include "LFDTD.h"


using namespace std;

const int nx(50), ny(50), nz(64), nx_TFSF{ 5 }, ny_TFSF{ 5 }, nz_TFSF{ 5 }, qstop(200), numProbe(1);
// const double t = 10e-9; // Simulation time duration
const double dt = 10e-12;
const int tStep = 1000;
int Nnode(0), NNZ(0);

double pi = M_PI;

double s = 6e10, fc(1e9), td(1 / (1.5 * fc)), tc(2e-9);
int pulseType = 3; // 1 - Gaussian, 2 - Gaussian Derivative, 3 - Modulated Gaussian

double eps0 = 8.854e-12;

double mu0 = 4 * pi * 1e-7;

double v0 = 1 / sqrt(eps0 * mu0);

double scale = 230 / (s * (tStep - 1) * dt / 2); // Similar as in Myunghyun's paper

// Current density source & resistor

double R0(50), Rs(0), Rl(0);
int jCellIndex[6] = { nx_TFSF + 1,nx - nx_TFSF,ny_TFSF + 1,ny + 1 - ny_TFSF,nz_TFSF + 1,nz + 1 - nz_TFSF }; // jCellIndex is assumed to be 1-indexed
int jCellIndex_Z[6] = { nx_TFSF + 1,nx + 1 - nx_TFSF,ny_TFSF + 1,ny + 1 - ny_TFSF,nz_TFSF + 1,nz - nz_TFSF };
int jDirecIndex[3] = { 1, 0, 1 };

int mCellIndex[6] = { nx_TFSF + 1,nx - nx_TFSF,ny_TFSF + 1,ny + 1 - ny_TFSF,nz_TFSF,nz - nz_TFSF + 1 }; // mCellIndex is assumed to be 1-indexed
int mCellIndex_Z[6] = { nx_TFSF + 1,nx - nx_TFSF,ny_TFSF,ny + 1 - ny_TFSF,nz_TFSF + 1,nz - nz_TFSF + 1 };
int mDirecIndex[3] = { 0, 1, 1 };

int rDirecIndex[3] = { 0, 0, 0 }; // Determine whether resistor/probe will be set-up or not (0 - no, 1 - yes) - Don't worry about this now
int jResistorIndex[6] = { nx / 2 + 15, nx / 2 + 15, ny / 2 + 1, ny / 2 + 1, 44, 44}; // jResistorIndex is assumed to be 1-indexed - Don't worry about this now

// Probe Definition

int probeCell[] = { nx / 2 + 15, nx / 2 + 15, ny / 2 + 1, ny / 2 + 1, 44, 44 }; //Don't worry about this now

int probeDirecIndex[3] = { 0,0,1 };// Don't worry about this now

int maxit = 1500;
double tol = 1e-6;

int main()

{
    Grid::Mesh_Grid("Mesh_MATLAB"); // Set the grid assignment for the structure

    LFDTD_Coe _LFDTD_coe;
    _LFDTD_coe.Coe_SET();

    LFDTD _LFDTD;
    try {
        _LFDTD.Solver_Select();
    }
    catch (int)
    {
        std::cerr << "Program EXIT!" << std::endl;
        exit(-1);
    }
    

    _LFDTD.PrintQ_set(1);
    _LFDTD.SparseA_COO(_LFDTD_coe);
    _LFDTD.COO2CSR();
    
    _LFDTD.Sim_start(_LFDTD_coe);
    
#if WRITE_PROBE == 1
    _LFDTD.result_write("TFSF_Probe_38_7.txt");
#endif

    return 0;
}

