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
#include <map>

using namespace std;

static const int N = 10000;
static char order[19] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
static char orderNew[19] = { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 };

std::vector<char> load_alphabet(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Unable to open alphabet file\n";
        return {};
    }

    std::vector<char> alphabet;
    char c;
    while (file.get(c)) {
        alphabet.push_back(c);
    }

    return alphabet;
}

std::vector<char> CreateFile(const std::vector<char>& alphabet, int process_id, int world_size) {
    srand(time(NULL) ^ process_id);

    int i;
    std::vector<int> entry(alphabet.size(), 0);
    std::vector<char> symbols;

    for (i = 0; i < N / world_size; ++i) {
        int index = rand() % alphabet.size();
        char symbol = alphabet[index];
        entry[index]++;
        symbols.push_back(symbol);
    }

    return symbols;
}

int UP(vector <double>& P, double q) {
// n-длинна обрабатываемой части массива
// q-вставляемая сумма
	int j = P.size(); //место вставляемого элемента
	vector<double> m_P;
	m_P.resize(j-1);
	for (int i = 0; i < P.size(); i++) {
		if (P[i] <= q) {
			j = i;
			break;
		}
	}
	for (int i = 0; i < j; i++) {
		m_P[i] = P[i];
	}
	m_P[j] = q;
	for (int i = j+1; i < P.size()-1; i++) {
		m_P[i] = P[i];
	}
	P = m_P;
	return j; 
}

void Down(vector<vector<int>>& C, vector<int>& len, int n, int j) {
// n - длинна обрабатываемой части массива
// j - номер буквы, которая временно исключается
	vector<int> c = C[j];
	int l = len[j];
	for (int i = j; i <= n - 2; i++) {
		if (C.size() <= i) C.push_back((vector<int>)0);
		C[i] = C[i+1];
		len[i] = len[i + 1];
	}
	C[n - 1] = c;
	C[n] = c;
	if (C[n - 1].size() <= l) C[n - 1].push_back(-1);
	if (C[n].size() <= l) C[n].push_back(-1);
	C[n - 1][l] = 0;
	C[n][l] = 1;
	len[n - 1] = l + 1;
	len[n] = l + 1;
}

void Huffman(vector<vector<int>>& C, vector<int>& len, vector<double>& P) {
	if (P.size() == 2) {
		C[0].push_back(0);
		C[1].push_back(1);
		len[0] = 1;
		len[1] = 1;
	}
	else {
		double q = P[P.size() - 2] + P[P.size()-1];
		int j = UP(P, q);
		int size_P = P.size();
		Huffman(C, len, P);
		Down(C, len, size_P, j);
	}
}

double Max(vector<double> p) {
	double max = 0;
	for (int i = 0; i < p.size(); i++) {
		if (p[i] >= max) max = p[i];
	}
	return max;
}

bool myfunction(double i, double j) { return (i > j); }

string CodingHuffman(string s_input, string s_output, vector<vector<int>> C) {
    string result;
    ifstream input(s_input);
    string string1; 
    ofstream output(s_output);
    
    while (getline(input, string1)) { 
        int n = 0;
        while (n != string1.size()) {
            for (int i = 0; i < C.size(); i++) {
                if (string1[n] == order[i]) {
                    for (int j = 0; j < C[i].size(); j++)
                        result += to_string(C[i][j]);  
                    n++;
                }
                
            }
        }
    }
    output << result;
    return result;
}

// Функция для подсчета вероятностей символов в файле
vector<double> compute_probabilities(const vector<char>& symbols) {
    map<char, int> counts;
    int total = 0;
    for (char c : symbols) {
        counts[c]++;
        total++;
    }
	
    vector<double> probabilities;
    for (auto& pair : counts) {
        probabilities.push_back(static_cast<double>(pair.second) / total);
    }
    return probabilities;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    std::vector<char> alphabet = load_alphabet("symbols.txt");
    if (alphabet.empty()) {
        MPI_Finalize();
        return 1;
    }
    std::vector<char> symbols = CreateFile(alphabet, world_rank, world_size);

    std::cout << "На узле " << world_rank << " сгенерировано " << symbols.size() << " символов.\n";

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

        std::cout << "Узел с rank 0 получил " << symbols.size() << " символов.\n";

        std::ofstream file("Library.txt");
        if (!file) {
            std::cerr << "Unable to open file for writing\n";
            MPI_Finalize();
            return 1;
        }

        for (char symbol : symbols) {
            file << symbol;
        }
        file.close();

        std::cout << "В файл записано " << symbols.size() << " символов.\n";

        vector<double> probabilities = compute_probabilities(symbols);
        vector<vector<int>> C(probabilities.size());
        vector<int> len(probabilities.size());
        Huffman(C, len, probabilities);
    }
    else {
        MPI_Send(symbols.data(), symbols.size(), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}

