#pragma once

#include "LFDTD_Coe.h"

// Define the format to printf MKL_INT values
#if !defined(MKL_ILP64)
#define IFORMAT "%i"
#else
#define IFORMAT "%lli"
#endif

class LFDTD
{
public:

    LFDTD(); //Constructor

    typedef void (LFDTD::*Sim_ptr)(LFDTD_Coe&); // PARDISO/CUDA solver ptr

    static void Solver_Select(); // Query the user to make a selection of CPU/GPU based solver

    void PrintQ_set(int); // Set whether to print out Q order on to screen

    void Sparse_A_Val(const std::string& path, int n, const int nnz); // This sets up the Sparse matrix entry value
    void SparseA_COO(const LFDTD_Coe& Coe);
    void COO2CSR();

    void Intel_PARDISO_TFSF(LFDTD_Coe& Coe);

    void Sim_start(LFDTD_Coe& Coe)
    {
        if (_Solver == _PARDISO)
        {
            sim_ptr = &LFDTD::Intel_PARDISO_TFSF;
        }

        (this->*sim_ptr)(Coe);
    }


    void result_write(const std::string& InputFile);
    void Eq_result_write(const std::string& InputFile);

    void result_write_app(const std::string& InputFile);
    void Eq_result_write_app(const std::string& InputFile);

    ~LFDTD()
    {} 

private:

    static Solver _Solver;

    std::vector<std::pair<int, int>> IA;
    std::vector<int> JA;
    std::vector<double> VAL;

    std::unique_ptr<double[]> a;
    std::unique_ptr<int[]> ja;
    std::unique_ptr<int[]> ia;

    std::unique_ptr<double[]> b;
    std::unique_ptr<double[]> bs;
    std::unique_ptr<double[]> x;

    std::unique_ptr<double[]> sumE;
    std::unique_ptr<double[]> lagPoly;
    std::unique_ptr<double[]> lagPoly_sum;
    std::unique_ptr<double[]> vtg;
    std::unique_ptr<double[]> probe;
    std::unique_ptr<double[]> recordEq;

    int PrintQ;
    int ex011, ex111, ey111, ez111, ex122, ex120, ex102, ex100, ex121, ex101, ex112, ex110, ex021, ex012;
    int ey101, ey212, ey012, ey210, ey010, ey112, ey110, ey201, ey211, ey011, ey102;
    int ez210, ez221, ez201, ez021, ez001, ez211, ez011, ez121, ez101, ez110, ez120;

    std::unique_ptr<int[]> ja_M;
    std::unique_ptr<int[]> ia_M;

    Sim_ptr sim_ptr;

};


