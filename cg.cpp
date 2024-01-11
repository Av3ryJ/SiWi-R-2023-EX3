#include <math.h>
#include <iostream>
#include <fstream>
#include "Timer.h"
#include <mpi.h>
#include <chrono>
#include <thread>

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
    for (int row = 1; row <= ny-1; row++) {
        for (int col = 1; col <= nx-1; col++) {
            result[row*rowlength+col] = f[row*rowlength + col]
                    -alpha*values[row*rowlength + col]
                    -gamma*values[(col-1) + row*rowlength]
                    -gamma*values[(col+1) + row*rowlength]
                    -beta*values[col + (row-1)*rowlength]
                    -beta*values[col + (row+1)*rowlength];
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
           z[row*rowlength+col]= alpha*values[row*rowlength + col]
                                    +gamma*values[(col-1) + row*rowlength]
                                    +gamma*values[(col+1) + row*rowlength]
                                    +beta*values[col + (row-1)*rowlength]
                                    +beta* values[col + (row+1)*rowlength];
        }
    }
}

void devide(int numberOfGridpoints, int pid, int N_P, int &first_row, int &number_of_rows) {
    number_of_rows = (numberOfGridpoints / N_P) + (pid == N_P-1 ? numberOfGridpoints % N_P : 0);
    first_row = pid * (int) (numberOfGridpoints / N_P);
}

int main(int argc, char* argv[]) {
    // MPI init and save own pid and number of processes for later
    MPI_Init(&argc, &argv);
    int total_number_of_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &total_number_of_processes);
    int pid;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

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
    int numberOfGridPoints = (nx+1)*(ny+1);
    //Vorfaktoren der Diskretisierung
    double alpha = 2/hx_squared + 2/hy_squared + 4*pi_squared;
    //Korrektur von A2 gamma & beta vertauscht
    double gamma = -1/hx_squared;
    double beta = -1/hy_squared;

    //Gitter mit Werten und Rand
    double *values = new double[numberOfGridPoints];
    initialize(nx, ny, values, hx);

    // richtige Funktionswerte
    double *f= new double[numberOfGridPoints];
    for(int y=0; y<=ny; y++){
        for(int x=0; x<=nx; x++){
            f[y*(nx+1)+x]= 4*pi_squared*sin(2*M_PI*x*hx)*sinh(2*M_PI*y*hy);
        }
    }
    // wait for all to finish startup
    MPI_Barrier(MPI_COMM_WORLD);

    int first;
    int len;
    devide(numberOfGridPoints, pid, total_number_of_processes, first, len);
    std::cout << "PID: " << pid << " first: " << first << " len: " << len << std::endl;

    //Timing start
    double time = 100.0;
    siwir::Timer timer;
    int cg_iterations = c;

    //Berechnung...
    double* residuum = new double[numberOfGridPoints];
    calculateResidualVector(values, f, nx, ny, alpha, beta, gamma, residuum);//r=f-A*u

    double delta0 = vectorDotProduct(residuum, residuum, numberOfGridPoints); //delt0= rt*r

    if (sqrt(delta0) > eps) { // Stop condition: ||r||<= eps
        //d=r
        double* d = new double[numberOfGridPoints];
        for (int i = 0; i < numberOfGridPoints; i++) {
            d[i] = residuum[i];
        }
        //TODO: residuum neu allokieren in "Streifen"-Größe
        for (int iteration = 0; iteration < c; iteration++) { //iterations

            double* z =  new double[numberOfGridPoints];
            // z = A*d
            stencilVectorMul(d, nx,ny, alpha, beta, gamma, z);
            // a = delt0/(dt*z)
            double a_zwischenergebnis = vectorDotProduct(d, z, numberOfGridPoints);
            //TODO: broadcast my 'zwischenergebnis' to all other processes
            //TODO: gather all zwischenergebnisse and add them up
            //      zwischenergebnis += irecive...
            double a = delta0 / a_zwischenergebnis;
            // u = u+a*d
            vectorPlusScaledVector(values, a, d, values, numberOfGridPoints);
            // r = r-a*z // each processes only calculates a part of r
            vectorPlusScaledVector(residuum, -a, z, residuum, numberOfGridPoints);
            //TODO: communicate u

            // delta1 = rt * r
            double delta1 = vectorDotProduct(residuum, residuum, numberOfGridPoints);
            //TODO: broadcast own delta1 sub-sum
            //TODO: gather delta1 sub-sums
            //      delta1 += irecv...

            // stop condition: ||r||<=eps
            if (sqrt(delta1) <= eps) {
                cg_iterations = iteration;
                delta0 = delta1;
                break;
            }
            // b = delta1/delta0
            double b = delta1/delta0;
            // d = r+b*d
            vectorPlusScaledVector(residuum, b, d, d, numberOfGridPoints);
            //TODO: d zusammenkleben msg_id = start index

            // delta0 = delta1
            delta0 = delta1;
        }
    }

    //Timing stoppen & ausgeben
    time = std::min(time, timer.elapsed());
    if(pid==0){
        std::cout << time << std::endl;
        std::cout << cg_iterations << std::endl;
        std::cout << sqrt(delta0) << std::endl;

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
            fileO << col*hx << " " << row*hy << " " << values[row*(nx+1) + col] << std::endl;
        }
    }
    fileO.close();
    }
    MPI_Finalize();
}