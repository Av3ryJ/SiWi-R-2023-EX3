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

double vectorDotProduct(double* vec1, double* vec2, int length, int start) {
    double sum = 0;
    for (int i = start; i < start+length; i++) {
        sum += vec1[i]*vec2[i];
    }
    return sum;
}

double allreduce_vectorDotProduct(double value_to_send) {
    double result;
    MPI_Allreduce(&value_to_send, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return result;
}

void vectorPlusScaledVector(double *vec1, double scalingFactor, double *vec2, double *outVec, int start_index, int length) {
    for (int i = start_index; i < start_index+length; i++) {
        outVec[i] = vec1[i] + scalingFactor*vec2[i];
    }
}

void stencilVectorMul(double* values, int nx, double alpha, double beta, double gamma, double *z, int start_index, int length){
    int rowlength = nx+1;
    for (int row = start_index; row < start_index+length; row++) {
        for (int col = 1; col < nx; col++) {
           z[row*rowlength+col]= alpha*values[row*rowlength + col]
                                    +gamma*values[(col-1) + row*rowlength]
                                    +gamma*values[(col+1) + row*rowlength]
                                    +beta*values[col + (row-1)*rowlength]
                                    +beta* values[col + (row+1)*rowlength];
        }
    }
}

void devide(int rows, int pid, int N_P, int &first_row, int &number_of_rows) {
    number_of_rows = (rows / N_P) + (pid == N_P - 1 ? rows % N_P : 0);
    first_row = pid * (int) (rows / N_P);
    if (pid == 0) {
        //to avoid boundary
        first_row++;
        number_of_rows--;
    }
    //to avoid boundary
    if (pid == N_P-1) number_of_rows--;
}

void devide_for_stitching(int rows, int pid, int N_P, int &first_row, int &number_of_rows) {
    number_of_rows = (rows / N_P) + (pid == N_P - 1 ? rows % N_P : 0);
    first_row = pid * (int) (rows / N_P);
    //now boundary is needed
}

void stitch_vector(double *vector, int N_P, int ny, int nx) {
    //Receive other vectors and stitch into vector
    for (int process = 0; process < N_P; process++) {
        //nothing to receive from self:
        int sender_len;
        int sender_start;
        devide_for_stitching(ny+1, process, N_P, sender_start, sender_len);
        MPI_Bcast(vector+sender_start*(nx+1), sender_len*(nx+1), MPI_DOUBLE, process, MPI_COMM_WORLD);
    }
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

    // nx und ny gibt Anzahl der Intervalle an, die entstehen → es gibt nx+1 & ny+1 Punkte in jede Richtung
    double hx = 2.0/nx; // length of one intervall in x-direction
    double hy = 1.0/ny; // length of one intervall in y-direction
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


    //Timing start
    double time = 100.0;
    siwir::Timer timer;
    int cg_iterations = c;

    //get process specific indices
    int first_index;
    int p_row_number;
    devide(ny+1, pid, total_number_of_processes, first_index, p_row_number);
    int len_p = p_row_number*(nx+1);

    //Berechnung...
    double* residuum = new double[numberOfGridPoints];
    calculateResidualVector(values, f, nx, ny, alpha, beta, gamma, residuum);//r=f-A*values

    double delta0 = vectorDotProduct(residuum, residuum, numberOfGridPoints, 0); //delta0= rt*r

    if (sqrt(delta0) > eps) { // Stop condition: ||r||<= eps
        //d=r
        double* d = new double[numberOfGridPoints];
        for (int i = 0; i < numberOfGridPoints; i++) {
            d[i] = residuum[i];
        }
        for (int iteration = 0; iteration < c; iteration++) { //iterations
            //each process only calculates a part of z
            // z = A*d
            double* z =  new double[numberOfGridPoints];
            stencilVectorMul(d, nx, alpha, beta, gamma, z,first_index, p_row_number);

            //only part of z -> part of a -> needs communication
            // a = delta0/(dt*z)
            double a_zwischenergebnis = vectorDotProduct(d, z, len_p, first_index*(nx+1));
            a_zwischenergebnis = allreduce_vectorDotProduct(a_zwischenergebnis);
            double a = delta0 / a_zwischenergebnis;

            // values = values+a*d
            vectorPlusScaledVector(values, a, d, values, 0, numberOfGridPoints);
            //TODO: communicate values only speedup not necessary

            // r = r-a*z // each processes only calculates a part of r
            vectorPlusScaledVector(residuum, -a, z, residuum, first_index*(nx+1), len_p);

            // delta1 = rt * r
            double delta1 = vectorDotProduct(residuum, residuum, len_p, first_index*(nx+1));
            //part of -> part of delta1 -> communication
            delta1 = allreduce_vectorDotProduct(delta1);

            // stop condition: ||r||<=eps
            if (sqrt(delta1) <= eps) {
                cg_iterations = iteration;
                delta0 = delta1;
                break;
            }
            // b = delta1/delta0
            double b = delta1/delta0;
            // d = r+b*d // only part from first_index to len_p is updated
            vectorPlusScaledVector(residuum, b, d, d, first_index*(nx+1), len_p);
            //d zusammenkleben
            stitch_vector(d, total_number_of_processes, ny, nx);
            // delta0 = delta1
            delta0 = delta1;
        }
    }

    //Timing stoppen & ausgeben
    time = std::min(time, timer.elapsed());
    if(pid==0){
        std::cout << time << std::endl;
        std::cout << cg_iterations << std::endl;
        //sqrt(delta0) = L2-Norm of residual
        std::cout << sqrt(delta0) << std::endl;

        // fuer gnuplot muss das so aussehen → Funktion schreiben, die solution.txt so schreibt
        // # x y u(x,y)
        //  0 0 0
        // .. .. ..

        //x und y Werte müssen in dem definierten Bereich liegen → [0,1] und[0,2]
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