#include <iostream>
#include <omp.h>
using namespace std;

// складываем векторы покоординатно
// глобальная переменная для всех потоков = shared
// локальная переменная, у каждого потока своя = private

int main() {
    int n = 70;
    int *a = new int[n];
    int *b = new int[n];
    int *c = new int[n];
    for(int i=0; i<n; i++){
        a[i] = 1;
        b[i] = 2;
        c[i] = 0;
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

    /*

#pragma omp parallel
    {
        int threads_count = omp_get_num_threads();
        int block_size = n / threads_count;
        for(int i = 0; i < threads_count; i++){
            for(int j = 0; j < block_size; j++) {
                c[block_size * i + j] = a[block_size * i + j] + b[block_size * i + j];
            }
            if (i == threads_count - 1){
                int j = block_size * (threads_count - 1);
                while(j < n){
                    c[j] = a[j] + b[j];
                    j++;
                }
            }
        }

    }

    cout << "c: ";
    for(int i=0; i<n; i++){
        cout << c[i] << " ";
    }
    cout << endl;
    */

    // omp сама производит распараллеливание
    // грубо говоря можно считать, что содержимое {} уже находится внутри цикла for по всем потокам

#pragma omp parallel
    {
        int threads_count = omp_get_num_threads();
        int thread_num = omp_get_thread_num();
        int block_size = n / threads_count;
        int shift = block_size * thread_num;
        for(int i = 0; i < block_size; i++) {
            c[shift + i] = a[shift + i] + b[shift + i];
        }

        if (thread_num == 0){
            int j = block_size * (threads_count - 1);
            while(j < n){
                c[j] = a[j] + b[j];
                j++;
            }
        }

    }

    cout << "c: ";
    for(int i=0; i<n; i++){
        cout << c[i] << " ";
    }
    cout << endl;

    return 0;
}
