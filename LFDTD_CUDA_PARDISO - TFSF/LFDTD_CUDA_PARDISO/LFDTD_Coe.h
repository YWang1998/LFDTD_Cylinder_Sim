#pragma once
#include "global.h"
#include "FileReader.h"

class Grid
{

public:

    Grid() {}; // Default constructor that does nothing
    static void Mesh_Grid(); // Allocate mesh profile for grid assignment
    static void Mesh_Grid(const std::string& InputFile); // Allocate mesh profile for grid assignment based on an input file
    ~Grid() // Destructor that delete all the dynamic allocated pointer variable
    {
        delete[] _dxh;
        delete[] _dyh;
        delete[] _dzh;
        delete[] _dxe;
        delete[] _dye;
        delete[] _dze;
    }
protected:

    static double_array_type _eps;
    static double_array_type _mu;

    static double_array_type _sigmax;
    static double_array_type _sigmay;
    static double_array_type _sigmaz;

    static double* _dxh;
    static double* _dyh;
    static double* _dzh;
    static double* _dxe;
    static double* _dye;
    static double* _dze;

};

class LFDTD_Coe: public Grid{

public:
    friend class LFDTD;
    LFDTD_Coe(); // Constructor
    void Coe_SET(); // Compute the Laguerre-FDTD coefficients
    ~LFDTD_Coe()
    {
        delete[] _waveform;
    }

private:
    node_array_type _nodeNum;

    double_array_type _cex;
    double_array_type _cey;
    double_array_type _cez;

    double_array_type _chx;
    double_array_type _chy;
    double_array_type _chz;

    double_array_type _Jx;
    double_array_type _Jy;
    double_array_type _Jz;
    double_array_type _Mx;
    double_array_type _My;
    double_array_type _Mz;
    double_array_type _Rz;

    double_array_type _hx;
    double_array_type _hy;
    double_array_type _hz;
    double_array_type _sumHx;
    double_array_type _sumHy;
    double_array_type _sumHz;

    double* _waveform;
};