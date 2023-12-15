#include <omp.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "Timer.h"


//Ausgabe der Matrix zur Überprüfung
void print_matrix(int nx, int ny, double *v) {
    for(int y=0; y<=ny; y++){
        for(int x=0; x<=nx; x++){
            std::cout << v[y*(nx+1)+x] << " ";
            }
        std::cout << std::endl;
    }
}

// Initialisieren der Matrix
void initialize(int nx, int ny, double *v, double hx){
    double sinhyb = sinh(2*M_PI);
    #pragma omp parallel for
    for(int y=0; y<=ny; ++y){
        for(int x=0; x<=nx; ++x){
            if(y==ny){
                v[y*(nx+1)+x]= sin(2*M_PI*x*hx)*sinhyb;
            }else{
                v[y*(nx+1)+x]= 0;
            }
        }
    }
}

// Berechnung der L2-Norm des Residuums
// Residuum = Abweichung der angenäherten Lösung von den "richtigen" Funktionswerten 
// Norm = Wurzel aus summierten, quadrierten Residuuen der einzelnen Punkte durch Anzahl aller inneren Punkte
double calculateResidual(double* values, double* f, int nx, int ny, double alpha, double beta, double gamma,int rowlength ) {
    double factor = 1.0/sqrt((nx-2)*(ny-2));
    double sum = 0;
    for (int row = 1; row < ny-1; row++) {
        for (int col = 1; col < nx-1; col++) {
            double residual = f[row*nx + col] - alpha*values[row*nx + col]
                            +gamma*values[(col-1) + row*rowlength]
                            +gamma*values[(col+1) + row*rowlength]
                            +beta*values[col + (row-1)*rowlength]
                            +beta* values[col + (row+1)*rowlength];
            sum += residual*residual;
        }
    }
    return factor*sqrt(sum);
}

int main(int argc, char* argv[]){
    if (argc < 5) {
        std::cout << "Usage: (mpirun -np <N>) ./cg <nx> <ny> <c> <eps>" << std::endl;
        return -1;
    }

    int nx = atoi(argv[1]);
    int ny = atoi(argv[2]);
    int c = atoi(argv[3]);
    double eps = atof(argv[4]);

    // nx und ny gibt Anzahl der Intervalle an, die entstehen -> es gibt nx+1 & ny+1 Punkte in jede Richtung
    double hx = 2.0/nx; // lenght of one intervall in x-direction
    double hy = 1.0/ny; // lenght of one intervall in y-direction
    double hx_squared = hx*hx;
    double hy_squared = hy*hy;
    double pi_squared = M_PI*M_PI;
    int rowlength = nx+1;

    //Timing start
    double time = 100.0;
    siwir::Timer timer;

    //Berechnung

    //Timing stoppen & ausgeben
    time = std::min(time, timer.elapsed());
    std::cout << time << std::endl;

    //Norm berechnen
    //double residual = calculateResidual(...);
    //std::cout << std::endl << "L2 Norm of the residual = " << residual << std::endl;

    // fuer gnuplot muss das so aussehen --> Funktion schreiben, die solution.txt so schreibt
    // # x y u(x,y)
    //  0 0 0 
    // .. .. .. 
    
    //x und y Werte müssen in dem definierten Bereich liegen -> [0,1] und[0,2]
    /*std::ofstream fileO ("solution.txt");
        fileO << "# x y u(x,y)"<< std::endl;
        for (int col = 0; col < nx+1; col++) {
            for (int row = 0; row < ny+1; row++) {
               fileO << col*hx << " " << row*hy << " " << values[row*(nx+1) + col] << std::endl;
            }
        }
        fileO.close();
    */
}