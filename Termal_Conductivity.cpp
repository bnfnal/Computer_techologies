#include <iostream>
#include <fstream>
#include <omp.h>

using namespace std;

// уравнение теплопроводности


int main() {

    setlocale(0, "");

    int N = 50;
    double T = 5;
    double Time = 200;
    double H = 1;
    double L = 1;
    double k = 1;
    double dx = H / N;
    double dy = L / N;
    double dt = (dx * dx) / 4;


    double **U_new = new double *[N];
    double **U_old = new double *[N];

    for (int i = 0; i < N; i += 1) {
        U_new[i] = new double[N];
        U_old[i] = new double[N];

        for (int j = 0; j < N; j += 1) {
            U_new[i][j] = 0;
            if (i == 0 || i == N - 1 || j == 0 || j == N - 1) {
                U_old[i][j] = T;
            } else {
                U_old[i][j] = 0;
            }
            cout << U_old[i][j] << " ";
        }
        cout << endl;
    }

    cout << endl << endl;

    int NT = (int) (Time / dt);
    double accuracy = 0.00001;
    int ind = 0;
    double max_U = 0;


    // последовательно

    time_t ts = time(NULL);

    for (int m = 0; m <= NT; m += 1) {
        ind = 0;
        max_U = 0;
        // parallel for + вывод в файл
        // как улучшить, чтобы не терять время
        for (int i = 0; i < N; i += 1) {
            for (int j = 0; j < N; j += 1) {
                if (i == 0 || i == N - 1 || j == 0 || j == N - 1) {
                    U_new[i][j] = U_old[i][j];
                } else {
                    U_new[i][j] = U_old[i][j] + dt * k * ((U_old[i + 1][j] - 2 * U_old[i][j] + U_old[i - 1][j]) / dx +
                                                          (U_old[i][j + 1] - 2 * U_old[i][j] + U_old[i][j - 1]) / dy);
//                    if (abs(U_new[i][j] - U_old[i][j]) / U_old[i][j] > max_U){
//                        max_U = abs(U_new[i][j] - U_old[i][j]) / U_old[i][j];
//                    }
                }
            }
        }

//        if (max_U < accuracy){
//            cout << endl << "accuracy = " << accuracy << ",  m = " << m << ",   Time = " << m*dt << ",   NT = " << NT << endl;
//            break;
//        }

        auto change = U_old;
        U_old = U_new;
        U_new = change;
    }

    ts = time(NULL) - ts;

    cout << "последовательно за время " << ts << endl;

    cout << endl;
    for (int i1 = 0; i1 < N; i1 += 1) {
        for (int j1 = 0; j1 < N; j1 += 1) {
            cout << U_old[i1][j1] << " ";
        }
        cout << endl;
    }
    cout << endl << endl;

    // параллельно

    // nowait

    for (int i = 0; i < N; i += 1) {
        for (int j = 0; j < N; j += 1) {
            U_new[i][j] = 0;
            if (i == 0 || i == N - 1 || j == 0 || j == N - 1) {
                U_old[i][j] = T;
            } else {
                U_old[i][j] = 0;
            }
        }
    }

    time_t tnw = time(NULL);

#pragma omp parallel
    for (int m = 0; m <= NT; m += 1) {
        ind = 0;
        max_U = 0;
#pragma omp for nowait
        for (int i = 0; i < N; i += 1) {
            for (int j = 0; j < N; j += 1) {
                if (i == 0 || i == N - 1 || j == 0 || j == N - 1) {
                    U_new[i][j] = U_old[i][j];
                } else {
                    U_new[i][j] = U_old[i][j] + dt * k * ((U_old[i + 1][j] - 2 * U_old[i][j] + U_old[i - 1][j]) / dx +
                                                          (U_old[i][j + 1] - 2 * U_old[i][j] + U_old[i][j - 1]) / dy);
//                  if (abs(U_new[i][j] - U_old[i][j]) / U_old[i][j] > max_U){
//                       max_U = abs(U_new[i][j] - U_old[i][j]) / U_old[i][j];
//                    }
                }
            }
        }

//        if (max_U < accuracy){
//            cout << endl << "accuracy = " << accuracy << ",  m = " << m << ",   Time = " << m*dt << ",   NT = " << NT << endl;
//            break;
//        }

        auto change = U_old;
        U_old = U_new;
        U_new = change;
    }

    tnw = time(NULL) - tnw;

    cout << "nowait за время " << tnw << endl;

    cout << endl;

    for(int i1 = 0; i1 < N; i1 += 1) {
        for (int j1 = 0; j1 < N; j1 += 1) {
            cout << U_old[i1][j1] << " ";
        }
        cout << endl;
    }


    // schedule(static)

    for(int i = 0; i < N; i += 1) {
        for(int j = 0; j < N; j += 1) {
            U_new[i][j] = 0;
            if (i == 0 || i == N-1 || j == 0 || j == N-1){
                U_old[i][j] = T;
            }
            else{
                U_old[i][j] = 0;
            }
        }
    }

    time_t tps = time(NULL);

    for(int m = 0; m <= NT; m += 1) {
        ind = 0;
        max_U = 0;
#pragma omp parallel for schedule(static)
        for(int i = 0; i < N; i += 1) {
            for(int j = 0; j < N; j += 1) {
                if (i == 0 || i == N-1 || j == 0 || j == N-1){
                    U_new[i][j] = U_old[i][j];
                }
                else{
                    U_new[i][j] = U_old[i][j] + dt*k*((U_old[i+1][j] - 2*U_old[i][j] + U_old[i-1][j])/dx + (U_old[i][j+1] - 2*U_old[i][j] + U_old[i][j-1])/dy);
//                    if (abs(U_new[i][j] - U_old[i][j]) / U_old[i][j] > max_U){
//                        max_U = abs(U_new[i][j] - U_old[i][j]) / U_old[i][j];
//                    }
                }
            }
        }

//        if (max_U < accuracy){
//            cout << endl << "accuracy = " << accuracy << ",  m = " << m << ",   Time = " << m*dt << ",   NT = " << NT << endl;
//            break;
//        }

        auto change = U_old;
        U_old = U_new;
        U_new = change;
    }

    tps = time(NULL) - tps;

    cout << "schedule(static) за время " << tps << endl;

    cout << endl;
    for(int i1 = 0; i1 < N; i1 += 1) {
        for (int j1 = 0; j1 < N; j1 += 1) {
            cout << U_old[i1][j1] << " ";
        }
        cout << endl;
    }

    // schedule(dynamic)

    for (int i = 0; i < N; i += 1) {
        for (int j = 0; j < N; j += 1) {
            U_new[i][j] = 0;
            if (i == 0 || i == N - 1 || j == 0 || j == N - 1) {
                U_old[i][j] = T;
            } else {
                U_old[i][j] = 0;
            }
        }
    }

    time_t tpd = time(NULL);

    for(int m = 0; m <= NT; m += 1) {
        ind = 0;
        max_U = 0;
#pragma omp parallel for schedule(dynamic)
        for(int i = 0; i < N; i += 1) {
            for(int j = 0; j < N; j += 1) {
                if (i == 0 || i == N-1 || j == 0 || j == N-1){
                    U_new[i][j] = U_old[i][j];
                }
                else{
                    U_new[i][j] = U_old[i][j] + dt*k*((U_old[i+1][j] - 2*U_old[i][j] + U_old[i-1][j])/dx + (U_old[i][j+1] - 2*U_old[i][j] + U_old[i][j-1])/dy);
//                    if (abs(U_new[i][j] - U_old[i][j]) / U_old[i][j] > max_U){
//                        max_U = abs(U_new[i][j] - U_old[i][j]) / U_old[i][j];
//                    }
                }
            }
        }

//        if (max_U < accuracy){
//            cout << endl << "accuracy = " << accuracy << ",  m = " << m << ",   Time = " << m*dt << ",   NT = " << NT << endl;
//            break;
//        }

        auto change = U_old;
        U_old = U_new;
        U_new = change;
    }

    tpd = time(NULL) - tpd;

    cout << "schedule(dynamic) за время " << tpd << endl;

    cout << endl;
    for(int i1 = 0; i1 < N; i1 += 1) {
        for (int j1 = 0; j1 < N; j1 += 1) {
            cout << U_old[i1][j1] << " ";
        }
        cout << endl;
    }

    cout << "schedule(static) за время " << tps << endl;

    cout << endl;
    for(int i1 = 0; i1 < N; i1 += 1) {
        for (int j1 = 0; j1 < N; j1 += 1) {
            cout << U_old[i1][j1] << " ";
        }
        cout << endl;
    }

    // shared(fronter)

    for (int i = 0; i < N; i += 1) {
        for (int j = 0; j < N; j += 1) {
            U_new[i][j] = 0;
            if (i == 0 || i == N - 1 || j == 0 || j == N - 1) {
                U_old[i][j] = T;
            } else {
                U_old[i][j] = 0;
            }
        }
    }

    time_t tsf = time(NULL);

    int fronter = 0;
    int i = 0;
    int treads_count = 0;

#pragma omp parallel private(i) shared(fronter)
    treads_count = omp_get_num_threads();
    fronter = treads_count;

    for(int m = 0; m <= NT; m += 1) {
        ind = 0;
        max_U = 0;
#pragma omp parallel private(i) shared(fronter)
        i = omp_get_thread_num();
        while(fronter <= treads_count) {
            for(int j = 0; j < N; j += 1) {
                if (i == 0 || i == N-1 || j == 0 || j == N-1){
                    U_new[i][j] = U_old[i][j];
                }
                else{
                    U_new[i][j] = U_old[i][j] + dt*k*((U_old[i+1][j] - 2*U_old[i][j] + U_old[i-1][j])/dx + (U_old[i][j+1] - 2*U_old[i][j] + U_old[i][j-1])/dy);
//                    if (abs(U_new[i][j] - U_old[i][j]) / U_old[i][j] > max_U){
//                        max_U = abs(U_new[i][j] - U_old[i][j]) / U_old[i][j];
//                    }
                }
            }
            i = fronter;
            fronter++;
        }

//        if (max_U < accuracy){
//            cout << endl << "accuracy = " << accuracy << ",  m = " << m << ",   Time = " << m*dt << ",   NT = " << NT << endl;
//            break;
//        }

        auto change = U_old;
        U_old = U_new;
        U_new = change;
    }

    tsf = time(NULL) - tsf;

    cout << "shared(fronter) за время " << tsf << endl;

    cout << endl;
    for(int i1 = 0; i1 < N; i1 += 1) {
        for (int j1 = 0; j1 < N; j1 += 1) {
            cout << U_old[i1][j1] << " ";
        }
        cout << endl;
    }

    cout << endl;
    cout << "последовательно за время " << ts << endl;
    cout << "nowait за время " << tnw << endl;
    cout << "schedule(static) за время " << tps << endl;
    cout << "schedule(dynamic) за время " << tpd << endl;
    cout << "shared(fronter) за время " << tsf << endl;

    return 0;
}

// чтобы ускорить процесс отключим синхронизацию по времени
// #pragma omp parallel for nowait
// но тогда потоки друг друга не ждут и работа заканчивается, когда закончил самый шустрый
// "гонка процессов"

// #pragma omp parallel for schedule(static)

// #pragma omp parallel for schedule(dynamic)
// если какой-то поток остановился, ему отдают задачи еще не закнчившего работу процесса
// это делается автоматически и занимает какое-то время
