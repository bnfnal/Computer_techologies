#include <iostream>
#include <omp.h>
using namespace std;

// складываем векторы покоординатно
// глобальная переменная для всех потоков = shared
// локальная переменная, у каждого потока своя = private

int main() {
    int n = 60;
    int tc = 8;
    int *a = new int[n];
    int *b = new int[n];
    for(int i=0; i<n; i++){
        a[i] = 1;
        b[i] = 2;
    }
    int c = 0;
    int d = 0;

    int *term = new int[tc];
    for(int i=0; i<tc; i++){
        term[i] = 0;
    }

    cout << "a: ";
    for(int i=0; i<n; i++){
        cout << a[i] << " ";
    }
    cout << endl;

    cout << "b: ";
    for(int i=0; i<n; i++){
        cout << b[i] << " ";
    }
    cout << endl;

// не работает, т.к. одновременно несколько потоков обращаеются к одной и той же переменной c

/*
 #pragma omp parallel
    {
        int threads_count = omp_get_num_threads();
        int thread_num = omp_get_thread_num();
        int block_size = n / threads_count;
        int shift = block_size * thread_num;

        for(int i = 0; i < block_size; i++) {
            c += a[shift + i] * b[shift + i];
        }

        if (thread_num == 0){
            int j = block_size * (threads_count - 1);
            while(j < n){
                c += a[j] * b[j];
                j++;
            }
        }

    }
*/
#pragma omp parallel
    {
        int threads_count = omp_get_num_threads();
        int thread_num = omp_get_thread_num();
        int block_size = n / threads_count;
        int shift = block_size * thread_num;

        for(int i = 0; i < block_size; i++) {
            term[thread_num] += a[shift + i] * b[shift + i];
        }

        if (thread_num == 0){
            int j = block_size * threads_count;
            while(j < n){
                term[0] += a[j] * b[j];
                j++;
            }
        }

    }

    for(int i = 0; i < tc; i++) {
        c += term[i];
    }

    cout << "c = " << c << endl;

    // встроенное распараллеливание

#pragma omp parallel for reduction (+:d)
    for(int i = 0; i < n; i++) {
        d += a[i] * b[i];
    }

    cout << "d = " << d << endl;

    return 0;
}
