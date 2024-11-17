//
// Created by Yifan Wang on 10/20/23.
//
#include "global.h"
#include "LFDTD_Coe.h"

double_array_type Grid::_eps;
double_array_type Grid::_mu;

double_array_type Grid::_sigmax;
double_array_type Grid::_sigmay;
double_array_type Grid::_sigmaz;

double* Grid::_dxh;
double* Grid::_dyh;
double* Grid::_dzh;
double* Grid::_dxe;
double* Grid::_dye;
double* Grid::_dze;


void Grid::Mesh_Grid(const std::string& InputFile)
{

    _eps.resize(boost::extents[nx + 1][ny + 1][nz + 1]);
    _mu.resize(boost::extents[nx + 1][ny + 1][nz + 1]);

    _sigmax.resize(boost::extents[nx][ny][nz]);
    _sigmay.resize(boost::extents[nx][ny][nz]);
    _sigmaz.resize(boost::extents[nx][ny][nz]);

    _dxh = new double[nx + 1];
    _dyh = new double[ny + 1];
    _dzh = new double[nz + 1];
    _dxe = new double[nx + 1];
    _dye = new double[ny + 1];
    _dze = new double[nz + 1];

    MatrixImport_1D(InputFile + "/dxh.txt", nx + 1, _dxh);
    MatrixImport_1D(InputFile + "/dyh.txt", ny + 1, _dyh);
    MatrixImport_1D(InputFile + "/dzh.txt", nz + 1, _dzh);

    MatrixImport_1D(InputFile + "/dxe.txt", nx + 1, _dxe);
    MatrixImport_1D(InputFile + "/dye.txt", ny + 1, _dye);
    MatrixImport_1D(InputFile + "/dze.txt", nz + 1, _dze);

    MatrixImport_3D(InputFile + "/sigmax_probe.txt", nx, ny, nz, _sigmax);
    MatrixImport_3D(InputFile + "/sigmay_probe.txt", nx, ny, nz, _sigmay);
    MatrixImport_3D(InputFile + "/sigmaz_probe.txt", nx, ny, nz, _sigmaz);

    fill_3D_array<double_array_type, double>(_eps, nx + 1, ny + 1, nz + 1, eps0);
    fill_3D_array<double_array_type, double>(_mu, nx + 1, ny + 1, nz + 1, mu0);
    
}

void Grid::Mesh_Grid()
{
    _eps.resize(boost::extents[nx + 1][ny + 1][nz + 1]);
    _mu.resize(boost::extents[nx + 1][ny + 1][nz + 1]);

    _sigmax.resize(boost::extents[nx][ny][nz]);
    _sigmay.resize(boost::extents[nx][ny][nz]);
    _sigmaz.resize(boost::extents[nx][ny][nz]);

    _dxh = new double[nx + 1];
    _dyh = new double[ny + 1];
    _dzh = new double[nz + 1];
    _dxe = new double[nx + 1];
    _dye = new double[ny + 1];
    _dze = new double[nz + 1];

    std::fill_n(_dxe, nx + 1, 0.009);
    std::fill_n(_dye, ny + 1, 0.009);
    std::fill_n(_dze, nz + 1, 0.01);

    std::fill_n(_dxh, nx + 1, 0.009);
    std::fill_n(_dyh, ny + 1, 0.009);
    std::fill_n(_dzh, nz + 1, 0.01);

    fill_3D_array<double_array_type, double>(_eps, nx + 1, ny + 1, nz + 1, eps0);
    fill_3D_array<double_array_type, double>(_mu, nx + 1, ny + 1, nz + 1, mu0);

}

LFDTD_Coe::LFDTD_Coe()
{

    _nodeNum.resize(boost::extents[nx + 1][ny + 1][nz + 1][3]);

    _cex.resize(boost::extents[nx + 1][ny + 1][nz + 1]);
    _cey.resize(boost::extents[nx + 1][ny + 1][nz + 1]);
    _cez.resize(boost::extents[nx + 1][ny + 1][nz + 1]);

    _chx.resize(boost::extents[nx + 1][ny + 1][nz + 1]);
    _chy.resize(boost::extents[nx + 1][ny + 1][nz + 1]);
    _chz.resize(boost::extents[nx + 1][ny + 1][nz + 1]);

    _waveform = new double[tStep];

    _Jx.resize(boost::extents[nx][ny + 1][nz + 1]);
    _Jy.resize(boost::extents[nx + 1][ny][nz + 1]);
    _Jz.resize(boost::extents[nx + 1][ny + 1][nz]);
    _Mx.resize(boost::extents[nx][ny][nz]);
    _My.resize(boost::extents[nx][ny][nz]);
    _Mz.resize(boost::extents[nx][ny][nz]);

    _Rz.resize(boost::extents[nx + 1][ny + 1][nz]);

    _hx.resize(boost::extents[nx][ny][nz]);
    _hy.resize(boost::extents[nx][ny][nz]);
    _hz.resize(boost::extents[nx][ny][nz]);
    _sumHx.resize(boost::extents[nx][ny][nz]);
    _sumHy.resize(boost::extents[nx][ny][nz]);
    _sumHz.resize(boost::extents[nx][ny][nz]);
}

void LFDTD_Coe::Coe_SET()
{

    for (int i = 0; i < nx + 1; ++i) {
        for (int j = 0; j < ny + 1; ++j) {
            for (int k = 0; k < nz + 1; ++k) {
                _cex[i][j][k] = 2 / (s * _eps[i][j][k] * _dxh[i]);
                _cey[i][j][k] = 2 / (s * _eps[i][j][k] * _dyh[j]);
                _cez[i][j][k] = 2 / (s * _eps[i][j][k] * _dzh[k]);
                _chx[i][j][k] = 2 / (s * _mu[i][j][k] * _dxe[i]);
                _chy[i][j][k] = 2 / (s * _mu[i][j][k] * _dye[j]);
                _chz[i][j][k] = 2 / (s * _mu[i][j][k] * _dze[k]);
            }
        }
    }

    Rl = R0 * (jResistorIndex[1] - jResistorIndex[0] + 1) * (jResistorIndex[3] - jResistorIndex[2] + 1) / (jResistorIndex[5] - jResistorIndex[4] + 1); // Load Resistance

    if (rDirecIndex[2] == 1)
    {
        for (int i = jResistorIndex[0] - 1; i <= jResistorIndex[1] - 1; ++i) {
            for (int j = jResistorIndex[2] - 1; j <= jResistorIndex[3] - 1; ++j) {
                for (int k = jResistorIndex[4] - 1; k <= jResistorIndex[5] - 1; ++k) {
                    _Rz[i][j][k] = 1;
                }
            }
        }

    }
    

    switch (pulseType) // Waveform Definition
    {
    case 1:
        for (int i = 0; i < tStep; ++i) {
            _waveform[i] = pow(exp(1), (-pow(((dt * (i + 1) - tc) / td), 2)));
        }
        break;

    case 3:
        for (int i = 0; i < tStep; ++i)
        {
            double coe = -(dt * (i + 1) - tc) * (dt * (i + 1) - tc) / (td * td);
            _waveform[i] = sin(2 * pi * fc * (dt * (i + 1) - tc)) * pow(exp(1), coe);
        }
        break;

    default:
        break;
    };

    for (int i = 0; i < nx + 1; ++i) {
        for (int j = 0; j < ny + 1; ++j) {
            for (int k = 0; k < nz + 1; ++k) {

                if (i != nx)
                {
                    _nodeNum[i][j][k][0] = Nnode;
                    Nnode++;
                }

                if (j != ny)
                {
                    _nodeNum[i][j][k][1] = Nnode;
                    Nnode++;
                }

                if (k != nz)
                {
                    _nodeNum[i][j][k][2] = Nnode;
                    Nnode++;
                }

            }
        }
    }

    printf("Laguerre-FDTD Meshing completed with total %d nodes.\n\n", Nnode);

}



