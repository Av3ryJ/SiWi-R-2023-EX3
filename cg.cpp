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

void calculateResidualVector(double *values, double *f, int nx, int ny, double alpha, double beta, double gamma, double *result) {
    int rowlength = nx+1;
    for (int row = 1; row < ny; row++) {
        for (int col = 1; col < nx; col++) {
            //warum result[row-1]?
            result[row-1] = f[row*nx + col] - alpha*values[row*nx + col]
                            +gamma*values[(col-1) + row*rowlength]
                            +gamma*values[(col+1) + row*rowlength]
                            +beta*values[col + (row-1)*rowlength]
                            +beta* values[col + (row+1)*rowlength];
        }
    }
}

double vectorDotProduct(double* vec1, double* vec2, int length) {
    double sum = 0;
    for (int i = 0; i < length; i++) {
        sum += vec1[i]*vec2[i];
    }
    return sum;
}

void vectorPlusScaledVector(double *vec1, double scalingFactor, double *vec2, double *outVec, int length) {
    for (int i = 0; i < length; i++) {
        outVec[i] = vec1[i] + scalingFactor*vec2[i];
    }
}

void stencilVectorMul(double* values, int nx, int ny, double alpha, double beta, double gamma, double *z){
    int rowlength = nx+1;
    
    for (int row = 1; row < ny; row++) {
        for (int col = 1; col < nx; col++) {
            //z[row-1] stimmt glaube nicht
           z[row-1]= alpha*values[row*nx + col]
                            +gamma*values[(col-1) + row*rowlength]
                            +gamma*values[(col+1) + row*rowlength]
                            +beta*values[col + (row-1)*rowlength]
                            +beta* values[col + (row+1)*rowlength];
        }
    }
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
    int numberOfInnerGridPoints = (nx-1)*(ny-1);
    //Vorfaktoren der Diskretisierung
    double alpha = 2/hx_squared + 2/hy_squared + 4*pi_squared;
    //Korrektur von A2 gamme & beta vertauscht
    double gamma = 1/hx_squared;
    double beta = 1/hy_squared;

    //Gitter mit Werten und Rand
    double *values = new double[(nx+1)*(ny+1)];
    initialize(nx, ny, values, hx);

    // richtige Funktionswerte
    double *f= new double[(nx+1)*(ny+1)];
    for(int y=0; y<=ny; y++){
        for(int x=0; x<=nx; x++){
            f[y*(nx+1)+x]= 4*pi_squared*sin(2*M_PI*x*hx)*sinh(2*M_PI*y*hy); 
        }
    }

    //Timing start
    double time = 100.0;
    siwir::Timer timer;

    //Berechnung...
    double* residuum = new double[numberOfInnerGridPoints]; 
    calculateResidualVector(values, f, nx, ny, alpha, beta, gamma, residuum);//r=f-A*u
    double delta0 = vectorDotProduct(residuum, residuum, numberOfInnerGridPoints); //delt0= rt*r
    if (sqrt(delta0) > eps) { // Stop condition: ||r||<= eps
        double* d = residuum; //d=r
        for (int count = 0; count < c; count++) { //iterations
            double* z =  new double[numberOfInnerGridPoints];
            stencilVectorMul(d, nx,ny, alpha, beta, gamma, z ); //z= A*d
            double a = delta0 / vectorDotProduct(d, z, numberOfInnerGridPoints); // a = delt0/(dt*z)
            vectorPlusScaledVector(values, a, d, values, numberOfInnerGridPoints); //RAND BEACHTEN!!!, u= u+a*d
            vectorPlusScaledVector(residuum, -a, z, residuum, numberOfInnerGridPoints);//r= r-a*z
            double delta1 = vectorDotProduct(residuum, residuum, numberOfInnerGridPoints);//delt1=rt*r
            if (sqrt(delta1) <= eps) break; //stop condition: ||r||<=eps
            double b = delta1/delta0; // b= delt1/delt2
            vectorPlusScaledVector(residuum, b, d, d, numberOfInnerGridPoints); //d= r+b*d
            delta0 = delta1; //delt0= delt1
        }
    }

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
    std::ofstream fileO ("solution.txt");
        fileO << "# x y u(x,y)"<< std::endl;
        for (int col = 0; col < nx+1; col++) {
            for (int row = 0; row < ny+1; row++) {
               //fileO << col*hx << " " << row*hy << " " << values[row*(nx+1) + col] << std::endl;
            }
        }
        fileO.close();
}