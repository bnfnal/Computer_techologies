#include <iostream>
#include <fstream>
#include <omp.h>

using namespace std;

// уравнение теплопроводности


int main() {
    int N = 10;
    double T = 5;
    double Time = 200;
    double H = 1;
    double L = 1;
    double k = 1;
    double dx = H/N;
    double dy = L/N;
    double dt = (dx * dx)/4;


    double **U_new = new double* [N];
    double **U_old = new double* [N];

    for(int i = 0; i < N; i += 1) {
        U_new[i] = new double[N];
        U_old[i] = new double[N];

        for(int j = 0; j < N; j += 1) {
            U_new[i][j] = 0;
            if (i == 0 || i == N-1 || j == 0 || j == N-1){
                U_old[i][j] = T;
            }
            else{
                U_old[i][j] = 0;
            }
            cout << U_old[i][j] << " ";
        }
        cout << endl;
    }

    int NT = (int)(Time/dt);
    double accuracy = 0.00001;
    int ind = 0;
    double max_U = 0;

    for(int m = 0; m <= NT; m += 1) {
        ind = 0;
        max_U = 0;
        // parallel for + вывод в файл
        // как улучшить, чтобы не терять время
        for(int i = 0; i < N; i += 1) {
            for(int j = 0; j < N; j += 1) {
                if (i == 0 || i == N-1 || j == 0 || j == N-1){
                    U_new[i][j] = U_old[i][j];
                }
                else{
                    U_new[i][j] = U_old[i][j] + dt*k*((U_old[i+1][j] - 2*U_old[i][j] + U_old[i-1][j])/dx + (U_old[i][j+1] - 2*U_old[i][j] + U_old[i][j-1])/dy);
                    if (abs(U_new[i][j] - U_old[i][j]) / U_old[i][j] > max_U){
                        max_U = abs(U_new[i][j] - U_old[i][j]) / U_old[i][j];
                    }
                }
            }
        }

        if (max_U < accuracy){
            cout << endl << "accuracy = " << accuracy << ",  m = " << m << ",   Time = " << m*dt << ",   NT = " << NT << endl;
            break;
        }

        auto change = U_old;
        U_old = U_new;
        U_new = change;
    }

    cout << endl;
    for(int i1 = 0; i1 < N; i1 += 1) {
        for (int j1 = 0; j1 < N; j1 += 1) {
            cout << U_old[i1][j1] << " ";
        }
        cout << endl;
    }

    return 0;
}
