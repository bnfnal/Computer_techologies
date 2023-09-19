#include <iostream>
#include <omp.h>
using namespace std;

// директива = #....
// любая директива omp.h = #pragma....
// то, что в {} выполняется параллельно, если нет {}, то параллельно выполняется то, что до первого ;
// нужна настройка компилятора
// в файл CMakeLists.txt нужно добавить строку set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fopenmp")
// кол-во надписей = кол-ву процессоров на компе (потоков)
// иногда надписи могут пресекаться и перемешиваться, это ок

// omp_get_num_threads();  // возвращ кол-во потоков
// omp_get_thread_num();  // возвращ номер активного потока (начало с 0)

int main() {
#pragma omp parallel
    {
        int thread_num = omp_get_num_threads();
        cout << "Hello, World! from " << omp_get_thread_num() << endl;
    }
    return 0;
}