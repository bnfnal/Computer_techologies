#include <iostream>
#include <fstream>
#include <omp.h>

using namespace std;

// уравнение теплопроводности


int main() {
    int N = 10;
    float T = 10;
    float H = 1;
    float L = 1;
    float k = 1;
    float dx = H/N;
    float dy = L/N;
    float dt = (dx * dx)/4;

    float **U_new = new float* [N];
    float **U_old = new float* [N];

    for(int i = 0; i < N; i += 1) {
        U_new[i] = new float[N];
        U_old[i] = new float[N];

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

    int Tn = (int)(T/dt);

    for(int m = 0; m <= 2; m += 1) {
        for(int i = 0; i < N; i += 1) {
            for(int j = 0; j < N; j += 1) {
                if (i == 0 || i == N-1 || j == 0 || j == N-1){
                    U_new[i][j] = U_old[i][j];
                }
                else{
                    U_new[i][j] = U_old[i][j] + dt*k*((U_old[i+1][j] - 2*U_old[i][j] + U_old[i-1][j])/dx + (U_old[i][j+1] - 2*U_old[i][j] + U_old[i][j-1])/dy);
                }
            }
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
