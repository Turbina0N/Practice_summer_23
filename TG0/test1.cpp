#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <string>
#include <algorithm>
#include <Windows.h>           
#include <sstream>

using namespace std;//предусмотрение возможности изменять размер массива arr
static const int N = 10000; //количество всех символов в файле
static char arr[] = { 'ж', 'з', 'и', 'к', 'л', 'м', ' ', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '@', '.' };
static char order[19] = { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 };
static char orderNew[19] = { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 };

// Генерация случайных символов. Она должна быть разной для каждого процесса,
// поэтому мы передаем номер процесса в функцию.
std::vector<char> CreateFile(int process_id, int world_size) {
    srand(process_id);  // Использование идентификатора процесса как базового значения.

    int i;
    std::vector<int> entry(19, 0);
    std::vector<double> probability(19, 0);
    std::vector<char> symbols;  // Возвращаемый вектор символов.

    for (i = 0; i < N / world_size; ++i)  // Генерируем N/world_size символов.
    {
        int index = rand() % 19;
        char symbol = arr[index];
        entry[index]++;
        symbols.push_back(symbol);  // Добавить символ в вектор
    }

    for (int i = 0; i < probability.size(); i++) {
        probability[i] = entry[i] / (double)(N / world_size);
    }

    return symbols;
}


int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_rank;  // номер этого процесса.
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
ъ
    int world_size;    //общее количество процессов.
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    std::vector<char> symbols = CreateFile(world_rank, world_size);
    if (world_rank == 0) {
        // Процесс 0 собирает все символы от всех процессов и записывает их в файл.
        for (int i = 1; i < world_size; ++i) {
            int count;
            MPI_Status status;

            // Сначала получить количество символов от другого процесса.
            MPI_Probe(i, 0, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_CHAR, &count);

            std::vector<char> other_symbols(count);
            MPI_Recv(other_symbols.data(), count, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //принимает сообщение с символами от процесса i и сохраняет его в вектор

            // Добавить полученные символы к нашим.
            symbols.insert(symbols.end(), other_symbols.begin(), other_symbols.end());
        }

        // Записать все символы в файл.
        std::ofstream file("symbols.txt");
        for (char symbol : symbols) {
            file << symbol;
        }
    }
    else {
        // Все остальные процессы отправляют свои символы процессу 0.
        MPI_Send(symbols.data(), symbols.size(), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}