#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>

using namespace std;

static const int N = 10000;
static char arr[] = { 'a', 'b', 'c', 'd', 'e', 'f', ' ', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '@', '.' };
static char order[19] = { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 };
static char orderNew[19] = { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 };

std::vector<char> CreateFile(int process_id, int world_size) {
    srand(process_id);

    int i;
    std::vector<int> entry(19, 0);
    std::vector<double> probability(19, 0);
    std::vector<char> symbols;

    for (i = 0; i < N / world_size; ++i) {
        int index = rand() % 19;
        char symbol = arr[index];
        entry[index]++;
        symbols.push_back(symbol);
    }

    for (int i = 0; i < probability.size(); i++) {
        probability[i] = entry[i] / (double)(N / world_size);
    }

    return symbols;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    std::vector<char> symbols = CreateFile(world_rank, world_size);
    if (world_rank == 0) {
        for (int i = 1; i < world_size; ++i) {
            int count;
            MPI_Status status;

            MPI_Probe(i, 0, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_CHAR, &count);

            std::vector<char> other_symbols(count);
            MPI_Recv(other_symbols.data(), count, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            symbols.insert(symbols.end(), other_symbols.begin(), other_symbols.end());
        }

        std::ofstream file("Library.txt");
        for (char symbol : symbols) {
            file << symbol;
        }
    }
    else {
        MPI_Send(symbols.data(), symbols.size(), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
