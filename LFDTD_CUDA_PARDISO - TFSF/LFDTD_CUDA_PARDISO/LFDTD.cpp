#include "global.h"
#include "LFDTD.h"

using namespace CBLAS;

Solver LFDTD::_Solver;

LFDTD::LFDTD()
{
    b = std::make_unique<double[]>(Nnode);
    bs = std::make_unique<double[]>(Nnode);
    x = std::make_unique<double[]>(Nnode);

    sumE = std::make_unique<double[]>(Nnode);
    lagPoly = std::make_unique<double[]>(4 * tStep);
    recordEq = std::make_unique<double[]>(numProbe * (qstop + 1));
    lagPoly_sum = std::make_unique<double[]>(tStep);
    vtg = std::make_unique<double[]>(tStep);
    probe = std::make_unique<double[]>(numProbe * tStep);

    PrintQ = 0;
}

void LFDTD::PrintQ_set(int i)
{
    if (i) PrintQ = 1;   
}

void LFDTD::Sparse_A_Val(const std::string& path, int n, const int nnz)
{

    ia = std::make_unique<int[]>(n+1);
    a = std::make_unique<double[]>(nnz);
    ja = std::make_unique<int[]>(nnz);

    Sparse_MatrixImport(path, n, nnz, a.get(), ia.get(), ja.get());
}

void LFDTD::SparseA_COO(const LFDTD_Coe& Coe)

{
    printf("Constructing Sparse Matrix A ...\n");
    // Ex equations
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 1; j < ny; ++j)
        {
            for (int k = 1; k < nz; ++k)
            {
                ex111 = Coe._nodeNum[i][j][k][0];

                ex121 = Coe._nodeNum[i][j + 1][k][0];
                ey211 = Coe._nodeNum[i + 1][j][k][1];
                ey111 = Coe._nodeNum[i][j][k][1];

                ex101 = Coe._nodeNum[i][j - 1][k][0];
                ey201 = Coe._nodeNum[i + 1][j - 1][k][1];
                ey101 = Coe._nodeNum[i][j - 1][k][1];

                ez211 = Coe._nodeNum[i + 1][j][k][2];
                ez111 = Coe._nodeNum[i][j][k][2];
                ex112 = Coe._nodeNum[i][j][k + 1][0];

                ez210 = Coe._nodeNum[i + 1][j][k - 1][2];
                ez110 = Coe._nodeNum[i][j][k - 1][2];
                ex110 = Coe._nodeNum[i][j][k - 1][0];

                IA.push_back({ NNZ,ex111 }); JA.emplace_back(ex111);
                VAL.emplace_back(1.0 + Coe._cey[i][j][k] * Coe._chy[i][j][k] + Coe._cey[i][j][k]
                    * Coe._chy[i][j - 1][k] + Coe._cez[i][j][k] * Coe._chz[i][j][k]
                    + Coe._cez[i][j][k] * Coe._chz[i][j][k - 1] + 2 * Coe._sigmax[i][j][k] / (s * Coe._eps[i][j][k]));
                NNZ++;

                IA.push_back({ NNZ,ex111 }); JA.emplace_back(ex121);
                VAL.emplace_back(-Coe._cey[i][j][k] * Coe._chy[i][j][k]); NNZ++;
                IA.push_back({ NNZ,ex111 }); JA.emplace_back(ey211);
                VAL.emplace_back(Coe._cey[i][j][k] * Coe._chx[i][j][k]); NNZ++;
                IA.push_back({ NNZ,ex111 }); JA.emplace_back(ey111);
                VAL.emplace_back(-Coe._cey[i][j][k] * Coe._chx[i][j][k]); NNZ++;

                IA.push_back({ NNZ,ex111 }); JA.emplace_back(ex101);
                VAL.emplace_back(-Coe._cey[i][j][k] * Coe._chy[i][j - 1][k]); NNZ++;
                IA.push_back({ NNZ,ex111 }); JA.emplace_back(ey201);
                VAL.emplace_back(-Coe._cey[i][j][k] * Coe._chx[i][j - 1][k]); NNZ++;
                IA.push_back({ NNZ,ex111 }); JA.emplace_back(ey101);
                VAL.emplace_back(Coe._cey[i][j][k] * Coe._chx[i][j - 1][k]); NNZ++;

                IA.push_back({ NNZ,ex111 }); JA.emplace_back(ez211);
                VAL.emplace_back(Coe._cez[i][j][k] * Coe._chx[i][j][k]); NNZ++;
                IA.push_back({ NNZ,ex111 }); JA.emplace_back(ez111);
                VAL.emplace_back(-Coe._cez[i][j][k] * Coe._chx[i][j][k]); NNZ++;
                IA.push_back({ NNZ,ex111 }); JA.emplace_back(ex112);
                VAL.emplace_back(-Coe._cez[i][j][k] * Coe._chz[i][j][k]); NNZ++;

                IA.push_back({ NNZ,ex111 }); JA.emplace_back(ez210);
                VAL.emplace_back(-Coe._cez[i][j][k] * Coe._chx[i][j][k - 1]); NNZ++;
                IA.push_back({ NNZ,ex111 }); JA.emplace_back(ez110);
                VAL.emplace_back(Coe._cez[i][j][k] * Coe._chx[i][j][k - 1]); NNZ++;
                IA.push_back({ NNZ,ex111 }); JA.emplace_back(ex110);
                VAL.emplace_back(-Coe._cez[i][j][k] * Coe._chz[i][j][k - 1]); NNZ++;

            }
        }
    }

    // Ey equations
    for (int i = 1; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int k = 1; k < nz; ++k)
            {
                ey111 = Coe._nodeNum[i][j][k][1];

                ey112 = Coe._nodeNum[i][j][k + 1][1];
                ez121 = Coe._nodeNum[i][j + 1][k][2];
                ez111 = Coe._nodeNum[i][j][k][2];

                ey110 = Coe._nodeNum[i][j][k - 1][1];
                ez120 = Coe._nodeNum[i][j + 1][k - 1][2];
                ez110 = Coe._nodeNum[i][j][k - 1][2];

                ex121 = Coe._nodeNum[i][j + 1][k][0];
                ex111 = Coe._nodeNum[i][j][k][0];
                ey211 = Coe._nodeNum[i + 1][j][k][1];

                ex021 = Coe._nodeNum[i - 1][j + 1][k][0];
                ex011 = Coe._nodeNum[i - 1][j][k][0];
                ey011 = Coe._nodeNum[i - 1][j][k][1];

                IA.push_back({ NNZ,ey111 }); JA.emplace_back(ey111);
                VAL.emplace_back(1.0 + Coe._cez[i][j][k] * Coe._chz[i][j][k] + Coe._cez[i][j][k]
                    * Coe._chz[i][j][k - 1] + Coe._cex[i][j][k] * Coe._chx[i][j][k]
                    + Coe._cex[i][j][k] * Coe._chx[i - 1][j][k] + 2 * Coe._sigmay[i][j][k] / (s * Coe._eps[i][j][k]));
                NNZ++;

                IA.push_back({ NNZ,ey111 }); JA.emplace_back(ey112);
                VAL.emplace_back(-Coe._cez[i][j][k] * Coe._chz[i][j][k]); NNZ++;
                IA.push_back({ NNZ,ey111 }); JA.emplace_back(ez121);;
                VAL.emplace_back(Coe._cez[i][j][k] * Coe._chy[i][j][k]); NNZ++;
                IA.push_back({ NNZ,ey111 }); JA.emplace_back(ez111);
                VAL.emplace_back(-Coe._cez[i][j][k] * Coe._chy[i][j][k]); NNZ++;

                IA.push_back({ NNZ,ey111 }); JA.emplace_back(ey110);
                VAL.emplace_back(-Coe._cez[i][j][k] * Coe._chz[i][j][k - 1]); NNZ++;
                IA.push_back({ NNZ,ey111 }); JA.emplace_back(ez120);
                VAL.emplace_back(-Coe._cez[i][j][k] * Coe._chy[i][j][k - 1]); NNZ++;
                IA.push_back({ NNZ,ey111 }); JA.emplace_back(ez110);
                VAL.emplace_back(Coe._cez[i][j][k] * Coe._chy[i][j][k - 1]); NNZ++;

                IA.push_back({ NNZ,ey111 }); JA.emplace_back(ex121);
                VAL.emplace_back(Coe._cex[i][j][k] * Coe._chy[i][j][k]); NNZ++;
                IA.push_back({ NNZ,ey111 }); JA.emplace_back(ex111);
                VAL.emplace_back(-Coe._cex[i][j][k] * Coe._chy[i][j][k]); NNZ++;
                IA.push_back({ NNZ,ey111 }); JA.emplace_back(ey211);
                VAL.emplace_back(-Coe._cex[i][j][k] * Coe._chx[i][j][k]); NNZ++;

                IA.push_back({ NNZ,ey111 }); JA.emplace_back(ex021);
                VAL.emplace_back(-Coe._cex[i][j][k] * Coe._chy[i - 1][j][k]); NNZ++;
                IA.push_back({ NNZ,ey111 }); JA.emplace_back(ex011);;
                VAL.emplace_back(Coe._cex[i][j][k] * Coe._chy[i - 1][j][k]); NNZ++;
                IA.push_back({ NNZ,ey111 }); JA.emplace_back(ey011);
                VAL.emplace_back(-Coe._cex[i][j][k] * Coe._chx[i - 1][j][k]); NNZ++;
            }
        }
    }

    // Ez equations
    for (int i = 1; i < nx; ++i)
    {
        for (int j = 1; j < ny; ++j)
        {
            for (int k = 0; k < nz; ++k)
            {
                if (Coe._Rz[i][j][k] == 1)
                {
                    ez111 = Coe._nodeNum[i][j][k][2];

                    ez211 = Coe._nodeNum[i + 1][j][k][2];
                    ex112 = Coe._nodeNum[i][j][k + 1][0];
                    ex111 = Coe._nodeNum[i][j][k][0];

                    ez011 = Coe._nodeNum[i - 1][j][k][2];
                    ex012 = Coe._nodeNum[i - 1][j][k + 1][0];
                    ex011 = Coe._nodeNum[i - 1][j][k][0];

                    ey112 = Coe._nodeNum[i][j][k + 1][1];
                    ey111 = Coe._nodeNum[i][j][k][1];
                    ez121 = Coe._nodeNum[i][j + 1][k][2];

                    ey102 = Coe._nodeNum[i][j - 1][k + 1][1];
                    ey101 = Coe._nodeNum[i][j - 1][k][1];
                    ez101 = Coe._nodeNum[i][j - 1][k][2];

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez111);;
                    VAL.emplace_back(1 + Coe._cex[i][j][k] * Coe._chx[i][j][k] + Coe._cex[i][j][k]
                        * Coe._chx[i - 1][j][k] + Coe._cey[i][j][k] * Coe._chy[i][j][k]
                        + Coe._cey[i][j][k] * Coe._chy[i][j - 1][k] + 2 * Coe._sigmaz[i][j][k] / (s * Coe._eps[i][j][k])
                        + 2 * Coe._dze[k] / (s * Coe._eps[i][j][k] * Coe._dxh[i] * Coe._dyh[j] * Rl));
                    NNZ++;

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez211);;
                    VAL.emplace_back(-Coe._cex[i][j][k] * Coe._chx[i][j][k]); NNZ++;
                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ex112);;
                    VAL.emplace_back(Coe._cex[i][j][k] * Coe._chz[i][j][k]); NNZ++;
                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ex111);;
                    VAL.emplace_back(-Coe._cex[i][j][k] * Coe._chz[i][j][k]); NNZ++;

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez011);;
                    VAL.emplace_back(-Coe._cex[i][j][k] * Coe._chx[i - 1][j][k]); NNZ++;
                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ex012);;
                    VAL.emplace_back(-Coe._cex[i][j][k] * Coe._chz[i - 1][j][k]); NNZ++;
                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ex011);;
                    VAL.emplace_back(Coe._cex[i][j][k] * Coe._chz[i - 1][j][k]); NNZ++;

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ey112);;
                    VAL.emplace_back(Coe._cey[i][j][k] * Coe._chz[i][j][k]); NNZ++;
                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ey111);;
                    VAL.emplace_back(-Coe._cey[i][j][k] * Coe._chz[i][j][k]); NNZ++;
                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez121);;
                    VAL.emplace_back(-Coe._cey[i][j][k] * Coe._chy[i][j][k]); NNZ++;

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ey102);;
                    VAL.emplace_back(-Coe._cey[i][j][k] * Coe._chz[i][j - 1][k]); NNZ++;
                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ey101);;
                    VAL.emplace_back(Coe._cey[i][j][k] * Coe._chz[i][j - 1][k]); NNZ++;
                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez101);;
                    VAL.emplace_back(-Coe._cey[i][j][k] * Coe._chy[i][j - 1][k]); NNZ++;

                }
                else
                {
                    ez111 = Coe._nodeNum[i][j][k][2];

                    ez211 = Coe._nodeNum[i + 1][j][k][2];
                    ex112 = Coe._nodeNum[i][j][k + 1][0];
                    ex111 = Coe._nodeNum[i][j][k][0];

                    ez011 = Coe._nodeNum[i - 1][j][k][2];
                    ex012 = Coe._nodeNum[i - 1][j][k + 1][0];
                    ex011 = Coe._nodeNum[i - 1][j][k][0];

                    ey112 = Coe._nodeNum[i][j][k + 1][1];
                    ey111 = Coe._nodeNum[i][j][k][1];
                    ez121 = Coe._nodeNum[i][j + 1][k][2];

                    ey102 = Coe._nodeNum[i][j - 1][k + 1][1];
                    ey101 = Coe._nodeNum[i][j - 1][k][1];
                    ez101 = Coe._nodeNum[i][j - 1][k][2];

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez111);;
                    VAL.emplace_back(1 + Coe._cex[i][j][k] * Coe._chx[i][j][k] + Coe._cex[i][j][k]
                        * Coe._chx[i - 1][j][k] + Coe._cey[i][j][k] * Coe._chy[i][j][k]
                        + Coe._cey[i][j][k] * Coe._chy[i][j - 1][k] + 2 * Coe._sigmaz[i][j][k] / (s * Coe._eps[i][j][k]));
                    NNZ++;

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez211);;
                    VAL.emplace_back(-Coe._cex[i][j][k] * Coe._chx[i][j][k]); NNZ++;
                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ex112);;
                    VAL.emplace_back(Coe._cex[i][j][k] * Coe._chz[i][j][k]); NNZ++;
                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ex111);;
                    VAL.emplace_back(-Coe._cex[i][j][k] * Coe._chz[i][j][k]); NNZ++;

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez011);;
                    VAL.emplace_back(-Coe._cex[i][j][k] * Coe._chx[i - 1][j][k]); NNZ++;
                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ex012);;
                    VAL.emplace_back(-Coe._cex[i][j][k] * Coe._chz[i - 1][j][k]); NNZ++;
                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ex011);;
                    VAL.emplace_back(Coe._cex[i][j][k] * Coe._chz[i - 1][j][k]); NNZ++;

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ey112);;
                    VAL.emplace_back(Coe._cey[i][j][k] * Coe._chz[i][j][k]); NNZ++;
                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ey111);;
                    VAL.emplace_back(-Coe._cey[i][j][k] * Coe._chz[i][j][k]); NNZ++;
                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez121);;
                    VAL.emplace_back(-Coe._cey[i][j][k] * Coe._chy[i][j][k]); NNZ++;

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ey102);;
                    VAL.emplace_back(-Coe._cey[i][j][k] * Coe._chz[i][j - 1][k]); NNZ++;
                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ey101);;
                    VAL.emplace_back(Coe._cey[i][j][k] * Coe._chz[i][j - 1][k]); NNZ++;
                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez101);;
                    VAL.emplace_back(-Coe._cey[i][j][k] * Coe._chy[i][j - 1][k]); NNZ++;
                }

            }
        }
    }

    // Outmost ABC boundary for Ex
    for (int i = 0; i < nx; ++i)
    {

        // Edge
        for (int j = 0; j < ny + 1; j += ny)
        {
            for (int k = 0; k < nz + 1; k += nz)
            {
                if ((j == 0) && (k == 0)) // Case (a-1)
                {
                    ex111 = Coe._nodeNum[i][j][k][0];
                    ex122 = Coe._nodeNum[i][j + 1][k + 1][0];

                    IA.push_back({ NNZ,ex111 }); JA.emplace_back(ex111);
                    VAL.emplace_back(1.0 / (sqrt(pow(Coe._dye[j], 2) + pow(Coe._dze[k], 2))) + s / (4.0 * v0)); NNZ++;

                    IA.push_back({ NNZ,ex111 }); JA.emplace_back(ex122);
                    VAL.emplace_back(-1.0 / (sqrt(pow(Coe._dye[j], 2) + pow(Coe._dze[k], 2))) + s / (4.0 * v0)); NNZ++;
                }
                else if ((j == 0) && (k == nz)) // Case (a-2)
                {
                    ex111 = Coe._nodeNum[i][j][k][0];
                    ex120 = Coe._nodeNum[i][j + 1][k - 1][0];

                    IA.push_back({ NNZ,ex111 }); JA.emplace_back(ex111);
                    VAL.emplace_back(1.0 / (sqrt(pow(Coe._dye[j], 2) + pow(Coe._dze[k - 1], 2))) + s / (4.0 * v0)); NNZ++;

                    IA.push_back({ NNZ,ex111 }); JA.emplace_back(ex120);
                    VAL.emplace_back(-1.0 / (sqrt(pow(Coe._dye[j], 2) + pow(Coe._dze[k - 1], 2))) + s / (4.0 * v0)); NNZ++;
                }
                else if ((j == ny) && (k == 0)) // Case (a-3)
                {
                    ex111 = Coe._nodeNum[i][j][k][0];
                    ex102 = Coe._nodeNum[i][j - 1][k + 1][0];

                    IA.push_back({ NNZ,ex111 }); JA.emplace_back(ex111);
                    VAL.emplace_back(1.0 / (sqrt(pow(Coe._dye[j - 1], 2) + pow(Coe._dze[k], 2))) + s / (4.0 * v0)); NNZ++;

                    IA.push_back({ NNZ,ex111 }); JA.emplace_back(ex102);
                    VAL.emplace_back(-1.0 / (sqrt(pow(Coe._dye[j - 1], 2) + pow(Coe._dze[k], 2))) + s / (4.0 * v0)); NNZ++;
                }
                else if ((j == ny) && (k == nz)) // Case (a-4)
                {
                    ex111 = Coe._nodeNum[i][j][k][0];
                    ex100 = Coe._nodeNum[i][j - 1][k - 1][0];

                    IA.push_back({ NNZ,ex111 }); JA.emplace_back(ex111);
                    VAL.emplace_back(1.0 / (sqrt(pow(Coe._dye[j - 1], 2) + pow(Coe._dze[k - 1], 2))) + s / (4.0 * v0)); NNZ++;

                    IA.push_back({ NNZ,ex111 }); JA.emplace_back(ex100);
                    VAL.emplace_back(-1.0 / (sqrt(pow(Coe._dye[j - 1], 2) + pow(Coe._dze[k - 1], 2))) + s / (4.0 * v0)); NNZ++;
                }
            }
        }

        // Face
        for (int j = 0; j < ny + 1; j += ny)
        {
            for (int k = 1; k < nz; ++k)
            {
                if (j == 0) // Case (1-1)
                {
                    ex111 = Coe._nodeNum[i][j][k][0];
                    ex121 = Coe._nodeNum[i][j + 1][k][0];

                    IA.push_back({ NNZ,ex111 }); JA.emplace_back(ex111);
                    VAL.emplace_back(1.0 / Coe._dye[j] + s / (4.0 * v0)); NNZ++;

                    IA.push_back({ NNZ,ex111 }); JA.emplace_back(ex121);
                    VAL.emplace_back(-1.0 / Coe._dye[j] + s / (4.0 * v0)); NNZ++;
                }
                else if (j == ny) // Case (1-2)
                {
                    ex111 = Coe._nodeNum[i][j][k][0];
                    ex101 = Coe._nodeNum[i][j - 1][k][0];

                    IA.push_back({ NNZ,ex111 }); JA.emplace_back(ex111);
                    VAL.emplace_back(1.0 / Coe._dye[j - 1] + s / (4.0 * v0)); NNZ++;

                    IA.push_back({ NNZ,ex111 }); JA.emplace_back(ex101);
                    VAL.emplace_back(-1.0 / Coe._dye[j - 1] + s / (4.0 * v0)); NNZ++;
                }
            }
        }

        for (int j = 1; j < ny; j++)
        {
            for (int k = 0; k < nz + 1; k += nz)
            {
                if (k == 0) // Case (1-3)
                {
                    ex111 = Coe._nodeNum[i][j][k][0];
                    ex112 = Coe._nodeNum[i][j][k + 1][0];

                    IA.push_back({ NNZ,ex111 }); JA.emplace_back(ex111);
                    VAL.emplace_back(1.0 / Coe._dze[k] + s / (4.0 * v0)); NNZ++;

                    IA.push_back({ NNZ,ex111 }); JA.emplace_back(ex112);
                    VAL.emplace_back(-1.0 / Coe._dze[k] + s / (4.0 * v0)); NNZ++;
                }
                else if (k == nz) // Case (1-4)
                {
                    ex111 = Coe._nodeNum[i][j][k][0];
                    ex110 = Coe._nodeNum[i][j][k - 1][0];

                    IA.push_back({ NNZ,ex111 }); JA.emplace_back(ex111);
                    VAL.emplace_back(1.0 / Coe._dze[k - 1] + s / (4.0 * v0)); NNZ++;

                    IA.push_back({ NNZ,ex111 }); JA.emplace_back(ex110);
                    VAL.emplace_back(-1.0 / Coe._dze[k - 1] + s / (4.0 * v0)); NNZ++;
                }
            }
        }
    }

    // Outmost ABC boundary for Ey
    for (int j = 0; j < ny; ++j)
    {

        // Edge
        for (int i = 0; i < nx + 1; i += nx)
        {
            for (int k = 0; k < nz + 1; k += nz)
            {
                if ((i == 0) && (k == 0)) // Case (b-1)
                {
                    ey111 = Coe._nodeNum[i][j][k][1];
                    ey212 = Coe._nodeNum[i + 1][j][k + 1][1];

                    IA.push_back({ NNZ,ey111 }); JA.emplace_back(ey111);
                    VAL.emplace_back(1.0 / (sqrt(pow(Coe._dxe[i], 2) + pow(Coe._dze[k], 2))) + s / (4.0 * v0)); NNZ++;

                    IA.push_back({ NNZ,ey111 }); JA.emplace_back(ey212);
                    VAL.emplace_back(-1.0 / (sqrt(pow(Coe._dxe[i], 2) + pow(Coe._dze[k], 2))) + s / (4.0 * v0)); NNZ++;
                }
                else if ((k == 0) && (i == nx)) // Case (b-2)
                {
                    ey111 = Coe._nodeNum[i][j][k][1];
                    ey012 = Coe._nodeNum[i - 1][j][k + 1][1];

                    IA.push_back({ NNZ,ey111 }); JA.emplace_back(ey111);
                    VAL.emplace_back(1.0 / (sqrt(pow(Coe._dxe[i - 1], 2) + pow(Coe._dze[k], 2))) + s / (4.0 * v0)); NNZ++;

                    IA.push_back({ NNZ,ey111 }); JA.emplace_back(ey012);
                    VAL.emplace_back(-1.0 / (sqrt(pow(Coe._dxe[i - 1], 2) + pow(Coe._dze[k], 2))) + s / (4.0 * v0)); NNZ++;
                }
                else if ((k == nz) && (i == 0)) // Case (b-3)
                {
                    ey111 = Coe._nodeNum[i][j][k][1];
                    ey210 = Coe._nodeNum[i + 1][j][k - 1][1];

                    IA.push_back({ NNZ,ey111 }); JA.emplace_back(ey111);
                    VAL.emplace_back(1.0 / (sqrt(pow(Coe._dxe[i], 2) + pow(Coe._dze[k - 1], 2))) + s / (4.0 * v0)); NNZ++;

                    IA.push_back({ NNZ,ey111 }); JA.emplace_back(ey210);
                    VAL.emplace_back(-1.0 / (sqrt(pow(Coe._dxe[i], 2) + pow(Coe._dze[k - 1], 2))) + s / (4.0 * v0)); NNZ++;
                }
                else if ((i == nx) && (k == nz)) // Case (b-4)
                {
                    ey111 = Coe._nodeNum[i][j][k][1];
                    ey010 = Coe._nodeNum[i - 1][j][k - 1][1];

                    IA.push_back({ NNZ,ey111 }); JA.emplace_back(ey111);
                    VAL.emplace_back(1.0 / (sqrt(pow(Coe._dxe[i - 1], 2) + pow(Coe._dze[k - 1], 2))) + s / (4.0 * v0)); NNZ++;

                    IA.push_back({ NNZ,ey111 }); JA.emplace_back(ey010);
                    VAL.emplace_back(-1.0 / (sqrt(pow(Coe._dxe[i - 1], 2) + pow(Coe._dze[k - 1], 2))) + s / (4.0 * v0)); NNZ++;
                }
            }
        }

        // Face
        for (int i = 1; i < nx; i++)
        {
            for (int k = 0; k < nz + 1; k += nz)
            {
                if (k == 0) // Case (2-1)
                {
                    ey111 = Coe._nodeNum[i][j][k][1];
                    ey112 = Coe._nodeNum[i][j][k + 1][1];

                    IA.push_back({ NNZ,ey111 }); JA.emplace_back(ey111);
                    VAL.emplace_back(1.0 / Coe._dze[k] + s / (4.0 * v0)); NNZ++;

                    IA.push_back({ NNZ,ey111 }); JA.emplace_back(ey112);
                    VAL.emplace_back(-1.0 / Coe._dze[k] + s / (4.0 * v0)); NNZ++;
                }
                else if (k == nz) // Case (2-2)
                {
                    ey111 = Coe._nodeNum[i][j][k][1];
                    ey110 = Coe._nodeNum[i][j][k - 1][1];

                    IA.push_back({ NNZ,ey111 }); JA.emplace_back(ey111);
                    VAL.emplace_back(1.0 / Coe._dze[k - 1] + s / (4.0 * v0)); NNZ++;

                    IA.push_back({ NNZ,ey111 }); JA.emplace_back(ey110);
                    VAL.emplace_back(-1.0 / Coe._dze[k - 1] + s / (4.0 * v0)); NNZ++;
                }
            }
        }

        for (int i = 0; i < nx + 1; i += nx)
        {
            for (int k = 1; k < nz; k++)
            {
                if (i == 0) // Case (2-3)
                {
                    ey111 = Coe._nodeNum[i][j][k][1];
                    ey211 = Coe._nodeNum[i + 1][j][k][1];

                    IA.push_back({ NNZ,ey111 }); JA.emplace_back(ey111);
                    VAL.emplace_back(1.0 / Coe._dxe[i] + s / (4.0 * v0)); NNZ++;

                    IA.push_back({ NNZ,ey111 }); JA.emplace_back(ey211);
                    VAL.emplace_back(-1.0 / Coe._dxe[i] + s / (4.0 * v0)); NNZ++;
                }
                else if (i == nx) // Case (2-4)
                {
                    ey111 = Coe._nodeNum[i][j][k][1];
                    ey011 = Coe._nodeNum[i - 1][j][k][1];

                    IA.push_back({ NNZ,ey111 }); JA.emplace_back(ey111);
                    VAL.emplace_back(1.0 / Coe._dxe[i - 1] + s / (4.0 * v0)); NNZ++;

                    IA.push_back({ NNZ,ey111 }); JA.emplace_back(ey011);
                    VAL.emplace_back(-1.0 / Coe._dxe[i - 1] + s / (4.0 * v0)); NNZ++;
                }
            }
        }
    }

    // Outmost ABC boundary for Ez
    for (int k = 0; k < nz; ++k) {

        // Edge

        for (int i = 0; i < nx + 1; i += nx) {
            for (int j = 0; j < ny + 1; j += ny) {
                if ((i == 0) && (j == 0)) // Case (c-1)
                {
                    ez111 = Coe._nodeNum[i][j][k][2];
                    ez221 = Coe._nodeNum[i + 1][j + 1][k][2];

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez111);
                    VAL.emplace_back(1.0 / (sqrt(pow(Coe._dxe[i], 2) + pow(Coe._dye[j], 2))) + s / (4.0 * v0)); NNZ++;

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez221);
                    VAL.emplace_back(-1.0 / (sqrt(pow(Coe._dxe[i], 2) + pow(Coe._dye[j], 2))) + s / (4.0 * v0)); NNZ++;

                }
                else if ((i == 0) && (j == ny)) // Case (c-2)
                {
                    ez111 = Coe._nodeNum[i][j][k][2];
                    ez201 = Coe._nodeNum[i + 1][j - 1][k][2];

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez111);
                    VAL.emplace_back(1.0 / (sqrt(pow(Coe._dxe[i], 2) + pow(Coe._dye[j - 1], 2))) + s / (4.0 * v0)); NNZ++;

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez201);
                    VAL.emplace_back(-1.0 / (sqrt(pow(Coe._dxe[i], 2) + pow(Coe._dye[j - 1], 2))) + s / (4.0 * v0)); NNZ++;

                }
                else if ((i == nx) && (j == 0)) // Case (c-3)
                {
                    ez111 = Coe._nodeNum[i][j][k][2];
                    ez021 = Coe._nodeNum[i - 1][j + 1][k][2];

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez111);
                    VAL.emplace_back(1.0 / (sqrt(pow(Coe._dxe[i - 1], 2) + pow(Coe._dye[j], 2))) + s / (4.0 * v0)); NNZ++;

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez021);
                    VAL.emplace_back(-1.0 / (sqrt(pow(Coe._dxe[i - 1], 2) + pow(Coe._dye[j], 2))) + s / (4.0 * v0)); NNZ++;

                }
                else if ((i == nx) && (j == ny)) // Case (c-4)
                {
                    ez111 = Coe._nodeNum[i][j][k][2];
                    ez001 = Coe._nodeNum[i - 1][j - 1][k][2];

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez111);
                    VAL.emplace_back(1.0 / (sqrt(pow(Coe._dxe[i - 1], 2) + pow(Coe._dye[j - 1], 2))) + s / (4.0 * v0)); NNZ++;

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez001);
                    VAL.emplace_back(-1.0 / (sqrt(pow(Coe._dxe[i - 1], 2) + pow(Coe._dye[j - 1], 2))) + s / (4.0 * v0)); NNZ++;

                }
            }
        }

        // Face
        for (int i = 0; i < nx + 1; i += nx) {
            for (int j = 1; j < ny; j++) {
                if (i == 0) // Case (3-1)
                {
                    ez111 = Coe._nodeNum[i][j][k][2];
                    ez211 = Coe._nodeNum[i + 1][j][k][2];

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez111);
                    VAL.emplace_back(1.0 / Coe._dxe[i] + s / (4.0 * v0)); NNZ++;

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez211);
                    VAL.emplace_back(-1.0 / Coe._dxe[i] + s / (4.0 * v0)); NNZ++;

                }
                else if (i == nx) // Case (3-2)
                {
                    ez111 = Coe._nodeNum[i][j][k][2];
                    ez011 = Coe._nodeNum[i - 1][j][k][2];

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez111);
                    VAL.emplace_back(1.0 / Coe._dxe[i - 1] + s / (4.0 * v0)); NNZ++;

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez011);
                    VAL.emplace_back(-1.0 / Coe._dxe[i - 1] + s / (4.0 * v0)); NNZ++;

                }
            }
        }

        for (int i = 1; i < nx; i++) {
            for (int j = 0; j < ny + 1; j += ny) {
                if (j == 0) // Case (3-3)
                {
                    ez111 = Coe._nodeNum[i][j][k][2];
                    ez121 = Coe._nodeNum[i][j + 1][k][2];

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez111);
                    VAL.emplace_back(1.0 / Coe._dye[j] + s / (4.0 * v0)); NNZ++;

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez121);
                    VAL.emplace_back(-1.0 / Coe._dye[j] + s / (4.0 * v0)); NNZ++;

                }
                else if (j == ny) // Case (3-4)
                {
                    ez111 = Coe._nodeNum[i][j][k][2];
                    ez101 = Coe._nodeNum[i][j - 1][k][2];

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez111);
                    VAL.emplace_back(1.0 / Coe._dye[j - 1] + s / (4.0 * v0)); NNZ++;

                    IA.push_back({ NNZ,ez111 }); JA.emplace_back(ez101);
                    VAL.emplace_back(-1.0 / Coe._dye[j - 1] + s / (4.0 * v0)); NNZ++;

                }
            }
        }
    }

}

void LFDTD::Solver_Select()
{

    int Solver_select;

    std::cout << "Please select the solver: [1] - PARDISO, [2] - cuBICGSTABL (Disabled) " << std::endl;
    std::cin >> Solver_select;

    if (Solver_select == 1)
    {
        LFDTD::_Solver = _PARDISO;
    }
    else
    {
        std::cerr << "No valid solver matches input selection!" << std::endl;;
        std::cerr << "Program exiting..." << std::endl;;
        throw -1;
    }

}

void LFDTD::COO2CSR()
{
    ia = std::make_unique<int[]>(Nnode + 1);
    a = std::make_unique<double[]>(NNZ);
    ja = std::make_unique<int[]>(NNZ);

    std::vector<std::pair<int, int>> JA_Group;
    std::vector<int> JA_Sorted_Idx;

    std::sort(IA.begin(), IA.end(), cmp);


    for (int i = 0; i < NNZ - 1; ++i)
    {

        if (IA[i].second == IA[i + 1].second)
        {
            JA_Group.push_back({ IA[i].first,JA[IA[i].first] });
        }
        else
        {
            JA_Group.push_back({ IA[i].first, JA[IA[i].first] });
            std::sort(JA_Group.begin(), JA_Group.end(), cmp);

            for (auto& ele : JA_Group)
            {
                JA_Sorted_Idx.push_back(ele.first);
            }

            JA_Group.clear();

        }

        ia[IA[i].second + 1]++;
    }

    // Account for the last element
    JA_Group.push_back({ IA[NNZ - 1].first, JA[IA[NNZ - 1].first] });
    std::sort(JA_Group.begin(), JA_Group.end(), cmp);

    for (auto& ele : JA_Group)
    {
        JA_Sorted_Idx.push_back(ele.first);
    }

    JA_Group.clear();
    ia[IA[NNZ - 1].second + 1]++;

    // 1-index in row column
    ia[0] = 1;

    for (int i = 0; i < Nnode; ++i)
    {
        ia[i + 1] += ia[i];
    }

    for (int i = 0; i < NNZ; ++i)
    {
        ja[i] = JA[JA_Sorted_Idx[i]] + 1;
        a[i] = VAL[JA_Sorted_Idx[i]];
    }
    

    IA.clear();
    JA.clear();
    VAL.clear();
    JA_Sorted_Idx.clear();

    printf("Sparse Matrix A Successfully Populated with CSR Format!\n");


}

void LFDTD::Intel_PARDISO_TFSF(LFDTD_Coe& Coe)
{
    int Pos, q(0);
    double jq(0);

    int nnode_1D = nz - nz_TFSF + 2;
    vec<int> nodeNum_1D = zeros<int>(nnode_1D);
    for (int i = 0; i < nnode_1D; ++i) nodeNum_1D[i] = i;

    int index_1D = 0;
    vec<int> row{ 13 * nnode_1D }; vec<int> col{ 13 * nnode_1D }; vec<double> val{ 13 * nnode_1D };

    for (int k = 1; k < nnode_1D - 1; ++k)
    {
        ex111 = nodeNum_1D[k];
        ex112 = nodeNum_1D[k+ 1];
        ex110 = nodeNum_1D[k- 1];

        row[index_1D] = ex111; col[index_1D] = ex111; val[index_1D] = std::round(1 + Coe._cez[0][0][k] * Coe._chz[0][0][k] + Coe._cez[0][0][k] * Coe._chz[0][0][k - 1]);
        ++index_1D;
        row[index_1D] = ex111; col[index_1D] = ex112; val[index_1D] = std::round(-Coe._cez[0][0][k] * Coe._chz[0][0][k]);
        ++index_1D;
        row[index_1D] = ex111; col[index_1D] = ex110; val[index_1D] = std::round(-Coe._cez[0][0][k] * Coe._chz[0][0][k - 1]);
        ++index_1D;
    }

    for (int k = nnode_1D - 1; k < nnode_1D; ++k)
    {
        ex111 = nodeNum_1D[k];
        ex110 = nodeNum_1D[k- 1];

        row[index_1D] = ex111; col[index_1D] = ex111; val[index_1D] = std::round(1 / Coe._dze[k - 1] + s / (4.0 * v0));
        ++index_1D;
        row[index_1D] = ex111; col[index_1D] = ex110; val[index_1D] = std::round(-1 / Coe._dze[k - 1] + s / (4.0 * v0));
        ++index_1D;
    }

    for (int k = 0; k < 1; ++k)
    {
        ex111 = nodeNum_1D[k];

        row[index_1D] = ex111; col[index_1D] = ex111; val[index_1D] = 1;
        ++index_1D;

    }

    row.resize(index_1D); col.resize(index_1D); val.resize(index_1D);
    mat<double> A_1D = LinMat(row, col, val);
    vec<double> x_1D{ nnode_1D }; vec<double> b_1D{ nnode_1D }; vec<double> sumE_1D{ nnode_1D };
    vec<double> hy_1D{ nnode_1D - 1 }; vec<double> sumHy_1D{ nnode_1D - 1 };

    auto start = std::chrono::high_resolution_clock::now();

    // PARDISO Set up
    MKL_INT mtype = 11;       /* Real unsymmetric matrix */
    // Descriptor of main sparse matrix properties
    struct matrix_descr descrA;
    // Structure with sparse matrix stored in CSR format
    sparse_matrix_t       csrA;
    /* RHS and solution vectors. */
    MKL_INT nrhs = 1;     /* Number of right hand sides. */
    /* Internal solver memory pointer pt, */
    /* 32-bit: int pt[64]; 64-bit: int pt[64] */
    /* or void *pt[64] should be OK on both architectures */
    void* pt[64];
    /* Pardiso control parameters. */
    MKL_INT iparm[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    /* Auxiliary variables. */
    double ddum;          /* Double dummy */
    MKL_INT idum;         /* Integer dummy. */
    /* -------------------------------------------------------------------- */
    /* .. Setup Pardiso control parameters. */
    /* -------------------------------------------------------------------- */
    for (int i = 0; i < 64; i++)
    {
        iparm[i] = 0;
    }
    iparm[0] = 1;         /* No solver default */
    iparm[1] = 2;         /* Fill-in reordering from METIS */
    iparm[3] = 0;         /* No iterative-direct algorithm */
    iparm[4] = 0;         /* No user fill-in reducing permutation */
    iparm[5] = 0;         /* Write solution into x */
    iparm[6] = 0;         /* Not in use */
    iparm[7] = 2;         /* Max numbers of iterative refinement stCoe._eps */
    iparm[8] = 0;         /* Not in use */
    iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;        /* Conjugate transposed/transpose solve */
    iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
    iparm[13] = 0;        /* Output: Number of perturbed pivots */
    iparm[14] = 0;        /* Not in use */
    iparm[15] = 0;        /* Not in use */
    iparm[16] = 0;        /* Not in use */
    iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1;       /* Output: Mflops for LU factorization */
    iparm[19] = 0;        /* Output: Numbers of CG Iterations */
    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;         /* Which factorization to use. */
    msglvl = 0;           /* Print statistical information  */
    error = 0;            /* Initialize error flag */
    /* -------------------------------------------------------------------- */
    /* .. Initialize the internal solver memory pointer. This is only */
    /* necessary for the FIRST call of the PARDISO solver. */
    /* -------------------------------------------------------------------- */
    for (int i = 0; i < 64; i++)
    {
        pt[i] = 0;
    }
    /* -------------------------------------------------------------------- */
    /* .. Reordering and Symbolic Factorization. This step also allocates */
    /* all memory that is necessary for the factorization. */
    /* -------------------------------------------------------------------- */
    phase = 11;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
        &Nnode, a.get(), ia.get(), ja.get(), &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

    if (error != 0)
    {
        printf("\nERROR during symbolic factorization: " IFORMAT, error);
        exit(1);
    }
    printf("\nReordering completed ... ");
    printf("\nNumber of nonzeros in factors = " IFORMAT, iparm[17]);
    printf("\nNumber of factorization MFLOPS = " IFORMAT, iparm[18]);
    /* -------------------------------------------------------------------- */
    /* .. Numerical factorization. */
    /* -------------------------------------------------------------------- */
    phase = 22;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
        &Nnode, a.get(), ia.get(), ja.get(), &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if (error != 0)
    {
        printf("\nERROR during numerical factorization: " IFORMAT, error);
        exit(2);
    }
    printf("\nFactorization completed ... \n");
    /* -------------------------------------------------------------------- */
    /* .. Back substitution and iterative refinement. */
    /* -------------------------------------------------------------------- */
    phase = 33;

    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
    descrA.mode = SPARSE_FILL_MODE_UPPER;
    descrA.diag = SPARSE_DIAG_NON_UNIT;
    mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ONE, Nnode, Nnode, ia.get(), ia.get() + 1, ja.get(), a.get());

    while (q <= qstop)
    {
        // Calculate Laguerre Polynomial
        if (q == 0)
        {
            for (int i = 0; i < tStep; ++i)
            {
                lagPoly[2 * tStep + i] = 1.0;
                lagPoly[3 * tStep + i] = -s * (i + 1) * dt / 2;
            }
        }
        else if (q == 1)
        {
            for (int i = 0; i < tStep; ++i)
            {
                lagPoly[1 * tStep + i] = lagPoly[2 * tStep + i];
                lagPoly[2 * tStep + i] = 1.0 - s * (i + 1) * dt;
            }
        }
        else
        {

            for (int i = 0; i < tStep; ++i)
            {
                lagPoly[0 * tStep + i] = lagPoly[1 * tStep + i];
                lagPoly[1 * tStep + i] = lagPoly[2 * tStep + i];
            }

            for (int i = 0; i < tStep; ++i)
            {
                lagPoly[2 * tStep + i] = (1.0 / q) * ((2.0 * q - 1.0 - s * (i + 1) * dt) * lagPoly[1 * tStep + i] - (q - 1) * lagPoly[0 * tStep + i]);
            }

            for (int i = 0; i < tStep; ++i)
            {
                if (lagPoly[2 * tStep + i] > 1e100) // Make sure that Laguerre polynomial does not go to infinity
                {
                    lagPoly[0 * tStep + i] = lagPoly[0 * tStep + i] * exp(-s * (i + 1) * dt / 2 * scale);
                    lagPoly[1 * tStep + i] = lagPoly[1 * tStep + i] * exp(-s * (i + 1) * dt / 2 * scale);
                    lagPoly[2 * tStep + i] = lagPoly[2 * tStep + i] * exp(-s * (i + 1) * dt / 2 * scale);
                    lagPoly[3 * tStep + i] = lagPoly[3 * tStep + i] + s * (i + 1) * dt / 2 * scale;
                }
            }

        }

        for (int i = 0; i < tStep; ++i)
        {
            lagPoly_sum[i] = lagPoly[2 * tStep + i] * exp(lagPoly[3 * tStep + i]);
        }

        // Compute Laguerre Coefficients for the source
        jq = 0;

        for (int i = 0; i < tStep; ++i)
        {
            jq += Coe._waveform[i] * (lagPoly[2 * tStep + i] * exp(lagPoly[3 * tStep + i])) * s * dt;
        }


        /* 1D TFSF Solving */

        /* 1D Ex equation */
        for (int k = 1; k < nnode_1D - 1; ++k)
        {
            ex111 = nodeNum_1D[k];

            b_1D[ex111] = 2 * Coe._cez[0][0][k] * (sumHy_1D[k] - sumHy_1D[k - 1]) - 2 * sumE_1D[ex111];
        }

        /* Outmost ABC boundary for Ex at z = Z */
        for (int k = nnode_1D - 1; k < nnode_1D; ++k)
        {
            ex111 = nodeNum_1D[k];
            ex110 = nodeNum_1D[k- 1];

            b_1D[ex111] = -s / (2 * v0) * (sumE_1D[ex111] + sumE_1D[ex110]);
        }

        /* Outmost TF/SF boundary */
        ex111 = nodeNum_1D[0];

        b_1D[ex111] = jq;

        /* Solve 1D x=A\b */

        A_1D.div(b_1D, x_1D);

        /* Update 1D variable */

        sumE_1D += x_1D;

        for (int k = 0; k < nnode_1D - 1; ++k)
        {
            ex111 = nodeNum_1D[k];
            ex112 = nodeNum_1D[k+ 1];

            hy_1D[k] = -Coe._chz[0][0][k] * (x_1D[ex112] - x_1D[ex111]) - 2 * sumHy_1D[k];
        }

        sumHy_1D += hy_1D;

        /* Current Source */
        
        // Jqx
        for (int i = jCellIndex[0] - 1; i < jCellIndex[1]; ++i)
        {
            for (int j = jCellIndex[2] - 1; j < jCellIndex[3]; j++)
            {
                // Bottom current source
                Coe._Jx[i][j][jCellIndex[4] - 1] = (-1/Coe._dzh[nz_TFSF-1])*hy_1D[0];

                // Top current source
                Coe._Jx[i][j][jCellIndex[5] - 1] = (1 / Coe._dzh[nz - nz_TFSF + 1]) * hy_1D[nz - 2*nz_TFSF + 1];
            }
        }

        // Left Jz current source
        for (int i = jCellIndex_Z[0] - 1; i < jCellIndex_Z[0]; ++i)
        {
            for (int j = jCellIndex_Z[2] - 1; j < jCellIndex_Z[3]; j++)
            {
                for (int k = jCellIndex_Z[4] - 1; k < jCellIndex_Z[5]; k++)
                {
                    Coe._Jz[i][j][k] = (1 / Coe._dxh[nz_TFSF - 1]) * hy_1D[k - nz_TFSF + 1];
                }
            }
        }

        // Right Jz current source
        for (int i = jCellIndex_Z[1] - 1; i < jCellIndex_Z[1]; ++i)
        {
            for (int j = jCellIndex_Z[2] - 1; j < jCellIndex_Z[3]; j++)
            {
                for (int k = jCellIndex_Z[4] - 1; k < jCellIndex_Z[5]; k++)
                {
                    Coe._Jz[i][j][k] = (-1 / Coe._dxh[nz_TFSF - 1]) * hy_1D[k - nz_TFSF + 1];
                }
            }
        }

        // Mqy
        for (int i = mCellIndex[0] - 1; i < mCellIndex[1]; ++i)
        {
            for (int j = mCellIndex[2] - 1; j < mCellIndex[3]; j++)
            {
                // Bottom magnetic source
                Coe._My[i][j][mCellIndex[4] - 1] = (-1 / Coe._dze[nz_TFSF - 1]) * x_1D[1];

                // Top magnetic source
                Coe._My[i][j][mCellIndex[5] - 1] = (1 / Coe._dze[nz - nz_TFSF + 1]) * x_1D[nnode_1D - nz_TFSF - 1];
            }
        }


        // Front magnetic source
        for (int i = mCellIndex_Z[0] - 1; i < mCellIndex_Z[1]; ++i)
        {
            for (int j = mCellIndex_Z[2] - 1; j < mCellIndex_Z[2]; j++)
            {
                for (int k = mCellIndex_Z[4] - 1; k < mCellIndex_Z[5]; k++)
                {
                    Coe._Mz[i][j][k] = (1 / Coe._dye[nz_TFSF - 1]) *x_1D[k - nz_TFSF + 1];
                }
            }
        }

        // Back magnetic source
        for (int i = mCellIndex_Z[0] - 1; i < mCellIndex_Z[1]; ++i)
        {
            for (int j = mCellIndex_Z[3] - 1; j < mCellIndex_Z[3]; j++)
            {
                for (int k = mCellIndex_Z[4] - 1; k < mCellIndex_Z[5]; k++)
                {
                    Coe._Mz[i][j][k] = (-1 / Coe._dye[nz_TFSF - 1]) * x_1D[k - nz_TFSF + 1];
                }
            }
        }

        // Build b vector

        // Ex equation except outmost PEC boundary
        // No re-assignment of b value for outmost PEC boundary

        for (int i = 0; i < nx; ++i)
        {
            for (int j = 1; j < ny; ++j)
            {
                for (int k = 1; k < nz; ++k)
                {
                    ex111 = Coe._nodeNum[i][j][k][0];
                    b[ex111] = -2.0 * Coe._cey[i][j][k] * (Coe._sumHz[i][j][k] - Coe._sumHz[i][j - 1][k]) + 2 * Coe._cez[i][j][k] * (Coe._sumHy[i][j][k] - Coe._sumHy[i][j][k - 1])
                        - 2.0 / (s * Coe._eps[i][j][k]) * Coe._Jx[i][j][k] + Coe._cez[i][j][k] * (2 / (s * Coe._mu[i][j][k])) * (Coe._My[i][j][k] - Coe._My[i][j][k-1]) \
                        - Coe._cey[i][j][k] * (2 / (s * Coe._mu[i][j][k])) * (Coe._Mz[i][j][k] - Coe._Mz[i][j-1][k]) \
                        - 2 * sumE[ex111];
                }
            }
        }

        // Ey equation except outmost PEC boundary
        // No re-assignment of b value for outmost PEC boundary

        for (int i = 1; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                for (int k = 1; k < nz; ++k)
                {
                    ey111 = Coe._nodeNum[i][j][k][1];
                    b[ey111] = -2.0 * Coe._cez[i][j][k] * (Coe._sumHx[i][j][k] - Coe._sumHx[i][j][k - 1]) + 2 * Coe._cex[i][j][k] * (Coe._sumHz[i][j][k] - Coe._sumHz[i - 1][j][k])
                        - 2.0 / (s * Coe._eps[i][j][k]) * Coe._Jy[i][j][k] - Coe._cez[i][j][k] * (2 / (s * Coe._mu[i][j][k])) * (Coe._Mx[i][j][k] - Coe._Mx[i][j][k-1]) \
                        + Coe._cex[i][j][k] * (2 / (s * Coe._mu[i][j][k])) * (Coe._Mz[i][j][k] - Coe._Mz[i-1][j][k])
                        - 2 * sumE[ey111];
                }
            }
        }

        // Ez equation except outmost PEC boundary
        // No re-assignment of b value for outmost PEC boundary

        for (int i = 1; i < nx; ++i)
        {
            for (int j = 1; j < ny; ++j)
            {
                for (int k = 0; k < nz; ++k)
                {
                    ez111 = Coe._nodeNum[i][j][k][2];
                    b[ez111] = -2.0 * Coe._cex[i][j][k] * (Coe._sumHy[i][j][k] - Coe._sumHy[i - 1][j][k])
                        + 2.0 * Coe._cey[i][j][k] * (Coe._sumHx[i][j][k] - Coe._sumHx[i][j - 1][k])
                        - 2.0 / (s * Coe._eps[i][j][k]) * Coe._Jz[i][j][k] + Coe._cey[i][j][k] * (2 / (s * Coe._mu[i][j][k])) * (Coe._Mx[i][j][k] - Coe._Mx[i][j-1][k]) \
                        - Coe._cex[i][j][k] * (2 / (s * Coe._mu[i][j][k])) * (Coe._My[i][j][k] - Coe._My[i-1][j][k])
                        - 2.0 * sumE[ez111];
                }
            }
        }

        // Outmost ABC boundary for Ex

        for (int i = 0; i < nx; ++i)
        {

            // Edge
            for (int j = 0; j < ny + 1; j += ny)
            {
                for (int k = 0; k < nz + 1; k += nz)
                {
                    if ((j == 0) && (k == 0)) // Case (a-1)
                    {
                        ex111 = Coe._nodeNum[i][j][k][0];
                        ex122 = Coe._nodeNum[i][j + 1][k + 1][0];
                        b[ex111] = -s / (2.0 * v0) * (sumE[ex111] + sumE[ex122]);
                    }
                    else if ((j == 0) && (k == nz)) // Case (a-2)
                    {
                        ex111 = Coe._nodeNum[i][j][k][0];
                        ex120 = Coe._nodeNum[i][j + 1][k - 1][0];
                        b[ex111] = -s / (2.0 * v0) * (sumE[ex111] + sumE[ex120]);
                    }
                    else if ((j == ny) && (k == 0)) // Case (a-3)
                    {
                        ex111 = Coe._nodeNum[i][j][k][0];
                        ex102 = Coe._nodeNum[i][j - 1][k + 1][0];
                        b[ex111] = -s / (2.0 * v0) * (sumE[ex111] + sumE[ex102]);
                    }
                    else if ((j == ny) && (k == nz)) // Case (a-4)
                    {
                        ex111 = Coe._nodeNum[i][j][k][0];
                        ex100 = Coe._nodeNum[i][j - 1][k - 1][0];
                        b[ex111] = -s / (2.0 * v0) * (sumE[ex111] + sumE[ex100]);
                    }
                }
            }

            // Face
            for (int j = 0; j < ny + 1; j += ny)
            {
                for (int k = 1; k < nz; ++k)
                {
                    if (j == 0) // Case (1-1)
                    {
                        ex111 = Coe._nodeNum[i][j][k][0];
                        ex121 = Coe._nodeNum[i][j + 1][k][0];
                        b[ex111] = -s / (2.0 * v0) * (sumE[ex111] + sumE[ex121]);
                    }
                    else if (j == ny) // Case (1-2)
                    {
                        ex111 = Coe._nodeNum[i][j][k][0];
                        ex101 = Coe._nodeNum[i][j - 1][k][0];
                        b[ex111] = -s / (2.0 * v0) * (sumE[ex111] + sumE[ex101]);
                    }
                }
            }

            for (int j = 1; j < ny; j++)
            {
                for (int k = 0; k < nz + 1; k += nz)
                {
                    if (k == 0) // Case (1-3)
                    {
                        ex111 = Coe._nodeNum[i][j][k][0];
                        ex112 = Coe._nodeNum[i][j][k + 1][0];
                        b[ex111] = -s / (2.0 * v0) * (sumE[ex111] + sumE[ex112]);
                    }
                    else if (k == nz) // Case (1-4)
                    {
                        ex111 = Coe._nodeNum[i][j][k][0];
                        ex110 = Coe._nodeNum[i][j][k - 1][0];
                        b[ex111] = -s / (2.0 * v0) * (sumE[ex111] + sumE[ex110]);
                    }
                }
            }
        }

        // Outmost ABC boundary for Ey

        for (int j = 0; j < ny; ++j)
        {

            // Edge
            for (int i = 0; i < nx + 1; i += nx)
            {
                for (int k = 0; k < nz + 1; k += nz)
                {
                    if ((i == 0) && (k == 0)) // Case (b-1)
                    {
                        ey111 = Coe._nodeNum[i][j][k][1];
                        ey212 = Coe._nodeNum[i + 1][j][k + 1][1];
                        b[ey111] = -s / (2.0 * v0) * (sumE[ey111] + sumE[ey212]);
                    }
                    else if ((k == 0) && (i == nx)) // Case (b-2)
                    {
                        ey111 = Coe._nodeNum[i][j][k][1];
                        ey012 = Coe._nodeNum[i - 1][j][k + 1][1];
                        b[ey111] = -s / (2.0 * v0) * (sumE[ey111] + sumE[ey012]);
                    }
                    else if ((k == nz) && (i == 0)) // Case (b-3)
                    {
                        ey111 = Coe._nodeNum[i][j][k][1];
                        ey210 = Coe._nodeNum[i + 1][j][k - 1][1];
                        b[ey111] = -s / (2.0 * v0) * (sumE[ey111] + sumE[ey210]);
                    }
                    else if ((i == nx) && (k == nz)) // Case (b-4)
                    {
                        ey111 = Coe._nodeNum[i][j][k][1];
                        ey010 = Coe._nodeNum[i - 1][j][k - 1][1];
                        b[ey111] = -s / (2.0 * v0) * (sumE[ey111] + sumE[ey010]);
                    }
                }
            }

            // Face
            for (int i = 1; i < nx; i++)
            {
                for (int k = 0; k < nz + 1; k += nz)
                {
                    if (k == 0) // Case (2-1)
                    {
                        ey111 = Coe._nodeNum[i][j][k][1];
                        ey112 = Coe._nodeNum[i][j][k + 1][1];
                        b[ey111] = -s / (2.0 * v0) * (sumE[ey111] + sumE[ey112]);
                    }
                    else if (k == nz) // Case (2-2)
                    {
                        ey111 = Coe._nodeNum[i][j][k][1];
                        ey110 = Coe._nodeNum[i][j][k - 1][1];
                        b[ey111] = -s / (2.0 * v0) * (sumE[ey111] + sumE[ey110]);
                    }
                }
            }

            for (int i = 0; i < nx + 1; i += nx)
            {
                for (int k = 1; k < nz; k++)
                {
                    if (i == 0) // Case (2-3)
                    {
                        ey111 = Coe._nodeNum[i][j][k][1];
                        ey211 = Coe._nodeNum[i + 1][j][k][1];
                        b[ey111] = -s / (2.0 * v0) * (sumE[ey111] + sumE[ey211]);
                    }
                    else if (i == nx) // Case (2-4)
                    {
                        ey111 = Coe._nodeNum[i][j][k][1];
                        ey011 = Coe._nodeNum[i - 1][j][k][1];
                        b[ey111] = -s / (2.0 * v0) * (sumE[ey111] + sumE[ey011]);
                    }
                }
            }
        }

        // Outmost ABC boundary for Ez

        for (int k = 0; k < nz; ++k)
        {

            // Edge

            for (int i = 0; i < nx + 1; i += nx)
            {
                for (int j = 0; j < ny + 1; j += ny)
                {
                    if ((i == 0) && (j == 0)) // Case (c-1)
                    {
                        ez111 = Coe._nodeNum[i][j][k][2];
                        ez221 = Coe._nodeNum[i + 1][j + 1][k][2];
                        b[ez111] = -s / (2.0 * v0) * (sumE[ez111] + sumE[ez221]);
                    }
                    else if ((i == 0) && (j == ny)) // Case (c-2)
                    {
                        ez111 = Coe._nodeNum[i][j][k][2];
                        ez201 = Coe._nodeNum[i + 1][j - 1][k][2];
                        b[ez111] = -s / (2.0 * v0) * (sumE[ez111] + sumE[ez201]);
                    }
                    else if ((i == nx) && (j == 0)) // Case (c-3)
                    {
                        ez111 = Coe._nodeNum[i][j][k][2];
                        ez021 = Coe._nodeNum[i - 1][j + 1][k][2];
                        b[ez111] = -s / (2.0 * v0) * (sumE[ez111] + sumE[ez021]);
                    }
                    else if ((i == nx) && (j == ny)) // Case (c-4)
                    {
                        ez111 = Coe._nodeNum[i][j][k][2];
                        ez001 = Coe._nodeNum[i - 1][j - 1][k][2];
                        b[ez111] = -s / (2.0 * v0) * (sumE[ez111] + sumE[ez001]);
                    }
                }
            }

            // Face
            for (int i = 0; i < nx + 1; i += nx)
            {
                for (int j = 1; j < ny; j++)
                {
                    if (i == 0) // Case (3-1)
                    {
                        ez111 = Coe._nodeNum[i][j][k][2];
                        ez211 = Coe._nodeNum[i + 1][j][k][2];
                        b[ez111] = -s / (2.0 * v0) * (sumE[ez111] + sumE[ez211]);
                    }
                    else if (i == nx) // Case (3-2)
                    {
                        ez111 = Coe._nodeNum[i][j][k][2];
                        ez011 = Coe._nodeNum[i - 1][j][k][2];
                        b[ez111] = -s / (2.0 * v0) * (sumE[ez111] + sumE[ez011]);
                    }
                }
            }

            for (int i = 1; i < nx; i++)
            {
                for (int j = 0; j < ny + 1; j += ny)
                {
                    if (j == 0) // Case (3-3)
                    {
                        ez111 = Coe._nodeNum[i][j][k][2];
                        ez121 = Coe._nodeNum[i][j + 1][k][2];
                        b[ez111] = -s / (2.0 * v0) * (sumE[ez111] + sumE[ez121]);
                    }
                    else if (j == ny) // Case (3-4)
                    {
                        ez111 = Coe._nodeNum[i][j][k][2];
                        ez101 = Coe._nodeNum[i][j - 1][k][2];
                        b[ez111] = -s / (2.0 * v0) * (sumE[ez111] + sumE[ez101]);
                    }
                }
            }
        }

        /* Call Intel PARDISO/CUDA BiCGSTABL solver here */

        // printf("\n\nSolving system with iparm[11] = " IFORMAT " ...\n", iparm[11]);
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
            &Nnode, a.get(), ia.get(), ja.get(), &idum, &nrhs, iparm, &msglvl, b.get(), x.get(), &error);
        if (error != 0)
        {
            printf("\nERROR during solution: " IFORMAT, error);
            exit(3);
        }

        // Compute residual

        /*
        mkl_sparse_d_mv(transA, 1.0, csrA, descrA, x, 0.0, bs);
        res = 0.0;
        res0 = 0.0;
        for (int j = 1; j <= Nnode; j++)
        {
            res += (bs[j - 1] - b[j - 1]) * (bs[j - 1] - b[j - 1]);
            res0 += b[j - 1] * b[j - 1];
        }
        res = sqrt(res) / sqrt(res0);
        // printf("\nRelative residual = %e", res);
        // Check residual
        if (res > 1e-10)
        {
            printf("Error: residual is too high!\n");
            exit(10);
        }
        */

        // Update sumE

        for (int i = 0; i < Nnode; ++i)
        {
            sumE[i] += x[i];
        }

        // Update Hx

        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                for (int k = 0; k < nz; ++k)
                {
                    ey112 = Coe._nodeNum[i][j][k + 1][1];
                    ey111 = Coe._nodeNum[i][j][k][1];
                    ez121 = Coe._nodeNum[i][j + 1][k][2];
                    ez111 = Coe._nodeNum[i][j][k][2];

                    Coe._hx[i][j][k] = Coe._chz[i][j][k] * (x[ey112] - x[ey111]) - Coe._chy[i][j][k] * (x[ez121] - x[ez111])
                        - 2.0 * Coe._sumHx[i][j][k] - (2.0 / (s * Coe._mu[i][j][k])) * Coe._Mx[i][j][k];
                    Coe._sumHx[i][j][k] += Coe._hx[i][j][k];
                }
            }
        }

        // Update Hy

        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                for (int k = 0; k < nz; ++k)
                {
                    ez211 = Coe._nodeNum[i + 1][j][k][2];
                    ez111 = Coe._nodeNum[i][j][k][2];
                    ex112 = Coe._nodeNum[i][j][k + 1][0];
                    ex111 = Coe._nodeNum[i][j][k][0];

                    Coe._hy[i][j][k] = Coe._chx[i][j][k] * (x[ez211] - x[ez111]) - Coe._chz[i][j][k] * (x[ex112] - x[ex111])
                        - 2.0 * Coe._sumHy[i][j][k] - (2.0 / (s * Coe._mu[i][j][k])) * Coe._My[i][j][k];
                    Coe._sumHy[i][j][k] += Coe._hy[i][j][k];
                }
            }
        }

        // Update Hz

        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                for (int k = 0; k < nz; ++k)
                {
                    ex121 = Coe._nodeNum[i][j + 1][k][0];
                    ex111 = Coe._nodeNum[i][j][k][0];
                    ey211 = Coe._nodeNum[i + 1][j][k][1];
                    ey111 = Coe._nodeNum[i][j][k][1];

                    Coe._hz[i][j][k] = Coe._chy[i][j][k] * (x[ex121] - x[ex111]) - Coe._chx[i][j][k] * (x[ey211] - x[ey111])
                        - 2.0 * Coe._sumHz[i][j][k] - (2.0 / (s * Coe._mu[i][j][k])) * Coe._Mz[i][j][k];
                    Coe._sumHz[i][j][k] += Coe._hz[i][j][k];
                }
            }
        }

        // Print the basis coefficient for the port with the lowest x and y index

        if (PrintQ == 1)
        {
            printf("q = %5d:", q);
            for (int i = 0; i < numProbe; ++i)
            {
                if (probeDirecIndex[0] == 1)
                {
                    Pos = Coe._nodeNum[probeCell[i * 6 + 0] - 1][probeCell[i * 6 + 2] - 1][probeCell[i * 6 + 4] - 1][0];
                    printf("p%d = %15.5e;", i + 1, x[Pos]);
                }
                else if (probeDirecIndex[1] == 1)
                {
                    Pos = Coe._nodeNum[probeCell[i * 6 + 0] - 1][probeCell[i * 6 + 2] - 1][probeCell[i * 6 + 4] - 1][1];
                    printf("p%d = %15.5e;", i + 1, x[Pos]);
                }
                else if (probeDirecIndex[2] == 1)
                {
                    Pos = Coe._nodeNum[probeCell[i * 6 + 0] - 1][probeCell[i * 6 + 2] - 1][probeCell[i * 6 + 4] - 1][2];
                    printf("p%d = %15.5e;", i + 1, x[Pos]);
                }

                else
                {
                    printf("Probe printing error\n");
                    exit(0);
                }
            }
            printf("\n");
        }

#if WRITE_PROBE

        for (int n = 0; n < numProbe; ++n)
        {
            std::fill_n(vtg.get(), tStep, 0.0);

            for (int i = probeCell[n * 6 + 0] - 1; i < probeCell[n * 6 + 1]; ++i)
            {
                for (int j = probeCell[n * 6 + 2] - 1; j < probeCell[n * 6 + 3]; ++j)
                {
                    for (int k = probeCell[n * 6 + 4] - 1; k < probeCell[n * 6 + 5]; ++k)
                    {
                        Pos = Coe._nodeNum[i][j][k][2];
                        for (int l = 0; l < tStep; ++l)
                        {
                            vtg[l] += x[Pos] * lagPoly_sum[l] * Coe._dze[k];
                            // probe[n][l] += vtg[l];
                        }
                    }
                }
            }
            for (int i = 0; i < tStep; i++)
            {
                probe[i + n * tStep] += vtg[i];
            }
            Pos = Coe._nodeNum[probeCell[n * 6 + 0] - 1][probeCell[n * 6 + 2] - 1][probeCell[n * 6 + 4] - 1][2];
            recordEq[q + n * (qstop + 1)] = x[Pos];
        }
#endif
        q++;
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::cout << "Elapsed time: " << elapsed_time << " milliseconds" << std::endl;

    /* -------------------------------------------------------------------- */
    /* .. Termination and release of memory. */
    /* -------------------------------------------------------------------- */
    phase = -1;           /* Release internal memory. */
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &Nnode, &ddum, ia.get(), ja.get(), &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

}

void LFDTD::result_write(const std::string& InputFile)
{
    std::ofstream myfile(InputFile);
    if (myfile.is_open())
    {
        for (int count = 0; count < tStep; count++)
        {
            for (int i = 0; i < numProbe; i++)
            {
                myfile << probe[count + i * tStep] << " ";
            }
            myfile << std::endl;
        }
        myfile.close();
    }

}

void LFDTD::Eq_result_write(const std::string& InputFile)
{
    std::ofstream myfile(InputFile);
    if (myfile.is_open())
    {
        for (int count = 0; count <= qstop; count++)
        {
            for (int i = 0; i < numProbe; i++)
            {
                myfile << recordEq[count + i * (qstop + 1)] << " ";
            }
            myfile << std::endl;
        }
        myfile.close();
    }
}

void LFDTD::result_write_app(const std::string& InputFile)
{
    std::ofstream myfile(InputFile, std::ios::app);
    if (myfile.is_open())
    {
        for (int count = 0; count < tStep; count++)
        {
            for (int i = 0; i < numProbe; i++)
            {
                myfile << probe[count + i * tStep] << " ";
            }
            myfile << std::endl;
        }
        myfile.close();
    }

    std::fill_n(probe.get(), tStep, 0.0);
}

void LFDTD::Eq_result_write_app(const std::string& InputFile)
{
    std::ofstream myfile(InputFile, std::ios::app);
    if (myfile.is_open())
    {
        for (int count = 0; count <= qstop; count++)
        {
            for (int i = 0; i < numProbe; i++)
            {
                myfile << recordEq[count + i * (qstop + 1)] << " ";
            }
            myfile << std::endl;
        }
        myfile.close();
    }
}