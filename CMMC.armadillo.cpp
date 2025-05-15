#include <iostream>
#include <armadillo>
#include <chrono>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
using namespace std;
using namespace arma;
/*
This code is the c++ version of Coupled matrix matrix completion algorithm. 
Armadillo(http://arma.sourceforge.net/) with OpenBLAS(https://www.openblas.net/) is needed to compile this code.
Complaing example:
g++ CMMC.cpp -o CMMC -std=c++11 -O3 -I/path/to/armadillo/<include>/ -I/path/to/openblas/<include>/ -L/path/to/openblas/openblas/<lib>/ -lopenblas
*/


vec generateOmegaAndV(mat &MXY, mat &Omega) {
    /*
    This function is to generate the Omega matrix and vector V.
    Omega is a two column matrix contains the row and column IDs for each known values in the main matrix.
    V is the vector contains all known values.
    */
    sp_mat mxy = sp_mat(MXY);
    vec ans(mxy.n_nonzero);
    Omega = mat(mxy.n_nonzero,2);
    
    typedef sp_mat::const_row_iterator iter;
    int k = 0;
    for (unsigned int i = 0; i < mxy.n_rows; i++) {
        for (iter j = mxy.begin_row(i); j != mxy.end_row(i); ++j) {
            ans.at(k++) = *j;
            //cout << i << "," << j.col() << endl;
            Omega((static_cast<__int64>(k)-1), 0) = int(i);
            Omega((static_cast<__int64>(k)-1), 1) = int(j.col());
        }
    }
    //cout << Omega.n_rows << "," << Omega.n_cols << endl;
    return ans;
}

void initializeZ(mat Omega, int m, cube &Z, field<mat> &B) {
    /* 
    This function is to initialize Z matrices. 
    Z matrices are position matrices that generated from Omega.
    For each coupled matrix, there will be a corresponding Z matrix.
    */
    int idx1, idy1 ,idx2, idy2; 
    for (int x = 0; x < m; ++x) {
        for (int y = 0; y < m; ++y) {
            idx1 = Omega(x, 0);
            idy1 = Omega(y, 0);
            idx2 = Omega(x, 1);
            idy2 = Omega(y, 1);
            Z(x, y, 0) = B(0, 0).at(idx1, idy1);
            Z(x, y, 1) = B(1, 0).at(idx2, idy2);
        }
    }
}

void updateZslice(vec Omega, int &m, mat &Z, mat &B) {
    /*
    In the main loop, Z matrices need to be updated for every iteration.
    This function is to update the Z matrices.
    */
    int idx, idy;
    for (int x = 0; x < m; ++x) {
        for (int y = 0; y < m; ++y) {
            idx = Omega(x);
            idy = Omega(y);
            Z(x, y) = B.at(idx,idy);
        }
    }
}

mat CoupledMatrixCompletion(mat &Omega, vec &v, field<mat> &C, int &iter_num) {
    /*
    This is the completion function.
    */
    // declare variables 
    int m = Omega.n_rows;
    int d = Omega.n_cols;
    int x, y;
    double dt;
    field<mat> B(2,1);
    cube Z(m, m, 2); // Z use cube
    mat b,M,E;
    mat Y(m,m);
    mat X(m,m);
    mat res; // final result matrix
    vec n(C.n_rows);
    
    Z.fill(0.0); // generate Z 
    //initialize B matrices, two B matrices will be initialized. 
    //these two matrices are used to generate the final result matirx.
    for (int i = 0; i < 2; ++i) {
        n(i) = C(i, 0).n_rows;
        B(i, 0) = b.eye(n(i), n(i));
    }
    
    initializeZ(Omega,m,Z,B); // filling values into Z based on Omega
    vec u = v;
    //main loop
    for (int l = 0; l < iter_num; ++l) { //main iteration, this is controled by the iter_num in main function
        for (int i = 0; i < 2; ++i) { //two coupled matrices will be used in here, there is a calculation for each one
            Y = Z.slice(!i); // select the corresponding Z matrix
            // generate M, M is the temporary matrix used to update B matrices.
            M.zeros(n(i), n(i));
            for (int j = 0; j < m; ++j) {
                for (int k = 0; k < m; ++k) {
                    x = Omega(j, i);
                    y = Omega(k, i);
                    M(x, y) = M(x, y) + u(j) * u(k) * Y(j,k);
                }
            }
            // update B
            B(i, 0) = C(i, 0) + B(i, 0) * M * B(i, 0);
            dt = pow(rcond(B(i, 0)),(-1/n(i)));
            B(i, 0) = dt * B(i, 0);
            // update Z
            updateZslice(Omega.col(i),m,Z.slice(i),B(i,0));
            X = Y % Z.slice(i);
            u = solve(X, v); // u is for final matrix
        }
    }
    // generate final result
    E.zeros(n(0),n(1));
    for (int i = 0; i < m; ++i) {
        E(Omega(i, 0), Omega(i, 1)) = u(i);
    }
    res = B(0, 0) * E * B(1,0);
    return(res);
}

int main(int argc, char** argv)
{
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    cout << "start" << endl;
    //declare vars
    // input files
    char *cgm_file;
    char *ggm_file;
    char *ccm_file;
    char *final_mat;

    ggm_file=argv[1];
    ccm_file=argv[2];
    cgm_file=argv[3];
    final_mat=argv[4];
    // declare mat objects
    mat MXX,MYY,MXY,Omega,W,D,I,MXY_hat;
    int a,b,s;
    // four hyperparameters
    double xscale = 0.9;
    double yscale = 0.9;
    double xepsilon = 0.9;
    double yepsilon = 0.9;

    int iter_num = 1; // interation number
    field<mat> C(2, 1); // declare a field var to hold two processed coupled matrices

    // load data from text file
    /*
    text file in tsv format without header or rownames.
    */
    MXX.load(ggm_file);
    MYY.load(ccm_file);
    MXY.load(cgm_file);

    W = MXY;
    a = MXX.n_rows;
    b = MYY.n_rows;
    s = Omega.n_rows;
    // generate Omega and V
    vec v = generateOmegaAndV(MXY, Omega);
    // preprocess two coupled matrices with four hyperparameters
    MXX = MXX * MXX;
    D = diagmat(pow(diagmat(MXX),-0.5));
    MXX = D * MXX * D;
    MYY = MYY * MYY;
    D = diagmat(pow(diagmat(MYY), -0.5));
    MYY = D * MYY * D;
    // use C to hold the processed coupled matrices
    C(0, 0) = xscale * MXX + xepsilon * I.eye(a, a);
    C(1, 0) = yscale * MYY + yepsilon * I.eye(b, b);
    // do the calculation
    MXY_hat = CoupledMatrixCompletion(Omega,v,C,iter_num);
    // save result into text file
    MXY_hat.save(final_mat,csv_ascii);

    cout << "Done" << endl;
    // print out the running time
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "Total elapsed time: " << chrono::duration_cast<chrono::seconds>(end - begin).count() << "[s]" << endl;
    return 0;
}

