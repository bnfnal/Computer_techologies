#include <iostream>
#include <omp.h>
using namespace std;

// последовательная работа потоков с помощью распараллеливания

int main() {
int i = 0;
#pragma omp parallel
    {
        int threads_count = omp_get_num_threads();
        while (i < threads_count) {
            if (omp_get_thread_num() == i) {
                cout << "Hello, World! from " << omp_get_thread_num() << endl;
                i++;
            }
        }
    }

    return 0;
}
