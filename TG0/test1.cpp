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

using namespace std;//ïðåäóñìîòðåíèå âîçìîæíîñòè èçìåíÿòü ðàçìåð ìàññèâà arr
static const int N = 10000; //êîëè÷åñòâî âñåõ ñèìâîëîâ â ôàéëå
static char arr[] = { 'a', 'b', 'c', 'd', 'e', 'f', ' ', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '@', '.' };
static char order[19] = { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 };
static char orderNew[19] = { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 };

// Ãåíåðàöèÿ ñëó÷àéíûõ ñèìâîëîâ. Îíà äîëæíà áûòü ðàçíîé äëÿ êàæäîãî ïðîöåññà,
// ïîýòîìó ìû ïåðåäàåì íîìåð ïðîöåññà â ôóíêöèþ.
std::vector<char> CreateFile(int process_id, int world_size) {
    srand(process_id);  // Èñïîëüçîâàíèå èäåíòèôèêàòîðà ïðîöåññà êàê áàçîâîãî çíà÷åíèÿ.

    int i;
    std::vector<int> entry(19, 0);
    std::vector<double> probability(19, 0);
    std::vector<char> symbols;  // Âîçâðàùàåìûé âåêòîð ñèìâîëîâ.

    for (i = 0; i < N / world_size; ++i)  // Ãåíåðèðóåì N/world_size ñèìâîëîâ.
    {
        int index = rand() % 19;
        char symbol = arr[index];
        entry[index]++;
        symbols.push_back(symbol);  // Äîáàâèòü ñèìâîë â âåêòîð
    }

    for (int i = 0; i < probability.size(); i++) {
        probability[i] = entry[i] / (double)(N / world_size);
    }

    return symbols;
}


int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_rank;  // íîìåð ýòîãî ïðîöåññà.
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
ú
    int world_size;    //îáùåå êîëè÷åñòâî ïðîöåññîâ.
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    std::vector<char> symbols = CreateFile(world_rank, world_size);
    if (world_rank == 0) {
        // Ïðîöåññ 0 ñîáèðàåò âñå ñèìâîëû îò âñåõ ïðîöåññîâ è çàïèñûâàåò èõ â ôàéë.
        for (int i = 1; i < world_size; ++i) {
            int count;
            MPI_Status status;

            // Ñíà÷àëà ïîëó÷èòü êîëè÷åñòâî ñèìâîëîâ îò äðóãîãî ïðîöåññà.
            MPI_Probe(i, 0, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_CHAR, &count);

            std::vector<char> other_symbols(count);
            MPI_Recv(other_symbols.data(), count, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //ïðèíèìàåò ñîîáùåíèå ñ ñèìâîëàìè îò ïðîöåññà i è ñîõðàíÿåò åãî â âåêòîð

            // Äîáàâèòü ïîëó÷åííûå ñèìâîëû ê íàøèì.
            symbols.insert(symbols.end(), other_symbols.begin(), other_symbols.end());
        }

        // Çàïèñàòü âñå ñèìâîëû â ôàéë.
        std::ofstream file("symbols.txt");
        for (char symbol : symbols) {
            file << symbol;
        }
    }
    else {
        // Âñå îñòàëüíûå ïðîöåññû îòïðàâëÿþò ñâîè ñèìâîëû ïðîöåññó 0.
        MPI_Send(symbols.data(), symbols.size(), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
