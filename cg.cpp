#include <cmath>
#include <iostream>
#include <fstream>
#include "Timer.h"

void print_matrix(int nx, int ny, double *v) {
    for (int y = 0; y <= ny; y++) {
        for (int x = 0; x <= nx; x++) {
            std::cout << v[y*(nx+1)+x] << " ";
        }
        std::cout << std::endl;
    }
}

void initialize(int nx, int ny, double *v, double hx) {
    double sinhyb = sinh(2 * M_PI);
    for (int y = 0; y <= ny; ++y) {
        for (int x = 0; x <= nx; ++x) {
            if (y == ny) {
                v[y*(nx+1)+x] = sin(2 * M_PI * x * hx) * sinhyb;
            } else {
                v[y*(nx+1)+x] = 0;
            }
        }
    }
}

void calculateResidualVector(double *values, double *f, int nx, int ny, double alpha, double beta, double gamma, double *result) {
    int rowlength = nx+1;
    for (int row = 1; row < ny; row++) {
        for (int col = 1; col < nx; col++) {
            int index = (row-1)*(nx-1) + (col-1);
            result[index] = f[row*rowlength + col] - alpha * values[row*rowlength + col]
                          + gamma * values[(col-1) + row*rowlength]
                          + gamma * values[(col+1) + row*rowlength]
                          + beta * values[col + (row-1)*rowlength]
                          + beta * values[col + (row+1)*rowlength];
        }
    }
}

double vectorDotProduct(double* vec1, double* vec2, int length) {
    double sum = 0;
    for (int i = 0; i < length; i++) {
        sum += vec1[i] * vec2[i];
    }
    return sum;
}

void vectorPlusScaledVector(double *vec1, double scalingFactor, double *vec2, double *outVec, int length) {
    for (int i = 0; i < length; i++) {
        outVec[i] = vec1[i] + scalingFactor * vec2[i];
    }
}

void stencilVectorMul(double* values, int nx, int ny, double alpha, double beta, double gamma, double *z) {
    int rowlength = nx+1;
    for (int row = 1; row < ny; row++) {
        for (int col = 1; col < nx; col++) {
            int index = (row-1)*(nx-1) + (col-1);
            z[index] = alpha * values[row*rowlength + col]
                     + gamma * values[(col-1) + row*rowlength]
                     + gamma * values[(col+1) + row*rowlength]
                     + beta * values[col + (row-1)*rowlength]
                     + beta * values[col + (row+1)*rowlength];
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cout << "Usage: ./cg <nx> <ny> <c> <eps>" << std::endl;
        return -1;
    }

    int nx = atoi(argv[1]);
    int ny = atoi(argv[2]);
    int c = atoi(argv[3]);
    double eps = atof(argv[4]);

    double hx = 2.0 / nx;
    double hy = 1.0 / ny;
    double pi_squared = M_PI * M_PI;
    int numberOfInnerGridPoints = (nx-1)*(ny-1);

    double alpha = 2 / (hx * hx) + 2 / (hy * hy) + 4 * pi_squared;
    double gamma = 1 / (hx * hx);
    double beta = 1 / (hy * hy);

    double *values = new double[(nx+1)*(ny+1)];
    double *f = new double[(nx+1)*(ny+1)];
    double *residuum = new double[numberOfInnerGridPoints];
    double *z = new double[numberOfInnerGridPoints];

    initialize(nx, ny, values, hx);
    for (int y = 0; y <= ny; y++) {
        for (int x = 0; x <= nx; x++) {
            f[y*(nx+1)+x] = 4 * pi_squared * sin(2 * M_PI * x * hx) * sinh(2 * M_PI * y * hy); 
        }
    }

    siwir::Timer timer;
    calculateResidualVector(values, f, nx, ny, alpha, beta, gamma, residuum);
    double delta0 = vectorDotProduct(residuum, residuum, numberOfInnerGridPoints);
    if (sqrt(delta0) > eps) {
        double* d = new double[numberOfInnerGridPoints];
        std::copy(residuum, residuum + numberOfInnerGridPoints, d);
        for (int count = 0; count < c; count++) {
            stencilVectorMul(d, nx, ny, alpha, beta, gamma, z);
            double a = delta0 / vectorDotProduct(d, z, numberOfInnerGridPoints);
            vectorPlusScaledVector(values, a, d, values, numberOfInnerGridPoints);
            vectorPlusScaledVector(residuum, -a, z, residuum, numberOfInnerGridPoints);
            double delta1 = vectorDotProduct(residuum, residuum, numberOfInnerGridPoints);
            if (sqrt(delta1) <= eps) break;
            double b = delta1 / delta0;
            vectorPlusScaledVector(residuum, b, d, d, numberOfInnerGridPoints);
            delta0 = delta1;
        }
        delete[] d;
    }

    double time = timer.elapsed();
    std::cout << time << std::endl;

    std::ofstream fileO("solution.txt");
    fileO << "# x y u(x,y)" << std::endl;
    for (int col = 0; col < nx+1; col++) {
        for (int row = 0; row < ny+1; row++) {
            fileO << col * hx << " " << row * hy << " " << values[row*(nx+1) + col] << std::endl;
        }
    }
    fileO.close();

    delete[] values;
    delete[] f;
    delete[] residuum;
    delete[] z;

    return 0;
}
