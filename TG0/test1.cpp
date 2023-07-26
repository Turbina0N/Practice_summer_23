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
static std::vector<char> order(19,-1);
static std::vector<char> orderNew(19,-1);

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
// std::cout << alphabet.size(); 
// for (auto c : alphabet) std::cout<<c<<" ";
    return alphabet;
}

std::vector<char> CreateFile(const std::vector<char>& alphabet, int process_id, int world_size) {
    srand(time(NULL) ^ process_id);
    
    int base_symbols_per_process = N / world_size;
    int remainder = N % world_size;
    int symbols_per_process = base_symbols_per_process + (process_id < remainder ? 1 : 0);

    std::vector<int> entry(alphabet.size(), 0);
    std::vector<char> symbols;

    for (int i = 0; i < symbols_per_process; ++i) {
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

// Функция для подсчета вероятностей символов в файле
vector<double> compute_probabilities(const vector<char>& symbols,const vector<char>& alphabet) {
    map<char, int> counts;
    int total = 0;
    for (char c : symbols) {
        counts[c]++;
        total++;
    }
    vector<double> probabilities;
    for (char c : alphabet) {
        if(counts.count(c) > 0) {
            probabilities.push_back(static_cast<double>(counts[c]) / total);
        } else {
            probabilities.push_back(0.0);
        }
    }
    return probabilities;
}
std::vector<std::vector<int>> transform_to_rectangle(const std::vector<std::vector<int>>& C) {
	int max_l=0;  
	for (size_t i = 0; i < C.size(); ++i) {
        if (C[i].size() > max_l) max_l = C[i].size();
    }
    std::vector<std::vector<int>> C_rectangular(C.size(), std::vector<int>(max_l, -1));
    for (size_t i = 0; i < C.size(); ++i) {
        copy(C[i].begin(), C[i].end(), C_rectangular[i].begin());
    }
    return C_rectangular;
}
// Вспомогательная функция для чтения файла в строку
std::string readFile(const std::string& fileName) {
    std::ifstream ifs(fileName.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
    std::ifstream::pos_type fileSize = ifs.tellg();
    ifs.seekg(0, std::ios::beg);

    std::vector<char> bytes(fileSize);
    ifs.read(&bytes[0], fileSize);

    return std::string(&bytes[0], fileSize);
}

std::string CodingHuffman(const std::string& s_input, const std::string& s_output, const std::vector<std::vector<int>>& C) {
    int rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Все процессы сначала читают весь файл
    std::string file_content = readFile(s_input);
    int total_symbols = file_content.size();

    int base_process = total_symbols / world_size;
    int remainder = total_symbols % world_size;

    // Определение начала и конца обработки каждого процесса
    int start_symbol = rank * base_process + std::min(rank, remainder);
    int symbols_per_process = base_process + (rank < remainder ? 1 : 0);
    int end_symbol = start_symbol + symbols_per_process;

    std::string result;
    for (int n = start_symbol; n < end_symbol; ++n) {
        for (int i = 0; i < C.size(); ++i) {
            if (file_content[n] == order[i]) {
                for (int j = 0; j < C[i].size(); ++j) {
                    if (C[i][j] != -1) {
                        result += std::to_string(C[i][j]);
                    }
                }
            }
        }
    }

    // // Запись результата в файл
    // std::ofstream output(s_output + "_part_" + std::to_string(rank)); // каждый процесс создает свой файл
    // output << result;

    return result;
}

std::string CodingRLE(std::string chunk) {
    std::string result;
    int count = 1;
    for (int i = 1; i < chunk.size(); i++) {
        if (chunk[i] == chunk[i - 1]) {
            count++;
        } else {
            if (chunk[i - 1] == '0') result += 'q';
            else if (chunk[i - 1] == '1') result += 'w';
            else if (chunk[i - 1] == '2') result += 'e';
            else if (chunk[i - 1] == '3') result += 'r';
            else if (chunk[i - 1] == '4') result += 't';
            else if (chunk[i - 1] == '5') result += 'y';
            else if (chunk[i - 1] == '6') result += 'u';
            else if (chunk[i - 1] == '7') result += 'i';
            else if (chunk[i - 1] == '8') result += 'o';
            else if (chunk[i - 1] == '9') result += 'p';
            else if (chunk[i - 1] == ' ') result += 'x';
            else result += chunk[i - 1];
            result += std::to_string(count);
            count = 1;
        }
        // Если текущий символ является последним в строке, то добавляем его и количество в результат
        if (i == chunk.size() - 1) {
            if (chunk[i] == '0') result += 'q';
            else if (chunk[i] == '1') result += 'w';
            else if (chunk[i] == '2') result += 'e';
            else if (chunk[i] == '3') result += 'r';
            else if (chunk[i] == '4') result += 't';
            else if (chunk[i] == '5') result += 'y';
            else if (chunk[i] == '6') result += 'u';
            else if (chunk[i] == '7') result += 'i';
            else if (chunk[i] == '8') result += 'o';
            else if (chunk[i] == '9') result += 'p';
            else if (chunk[i] == ' ') result += 'x';
            else result += chunk[i];
            result += std::to_string(count);
        }
    }
    return result;
}

void CodingRLE_MPI(string s_input, string s_output) {
    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    string file_content = readFile(s_input);

    int total_symbols = file_content.size();
    int base_process = total_symbols / world_size;
    int remainder = total_symbols % world_size;
    int start_symbol = world_rank * base_process + std::min(world_rank, remainder);
    int symbols_per_process = base_process + (world_rank < remainder ? 1 : 0);

    string chunk = file_content.substr(start_symbol, symbols_per_process);
    string result = CodingRLE(chunk);

    if (world_rank == 0) {
        ofstream output(s_output);
        output << result;
        for (int i = 1; i < world_size; i++) {
            int recv_size;
            MPI_Recv(&recv_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            char* recv_result = new char[recv_size];
            MPI_Recv(recv_result, recv_size, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            output << recv_result;
            delete[] recv_result;
        }
        output.close();
    } else {
        int send_size = result.size() + 1;
        MPI_Send(&send_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        const char *send_result = result.c_str();
        MPI_Send(send_result, send_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
}














int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    //std::vector<char> alphabet = load_alphabet("symbols.txt");
    std::vector<char> alphabet = { 'q', 'w', 'e', 'r', 't', 'y', ' ', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '@', '.' };
    std::cout <<alphabet.size();
    if (alphabet.empty()) {
        MPI_Finalize();
        return 1;
    }

    std::vector<char> symbols = CreateFile(alphabet, world_rank, world_size);
	
 //    std::cout <<"rank = "<< world_rank << "\t"<<std::endl;
 //    for (auto c: symbols){
	//     std::cout << c;
 //    }
	// std::cout<<std::endl;
	
    int numRows = 0;
    int numCols = 0;
    std::vector<std::vector<int>> C_rectangular;
	
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
	// конец генерации !!


	    
 	std::vector<double> probabilities = compute_probabilities(symbols, alphabet);
        std::vector<vector<int>> C(probabilities.size());
	std::cout << C.size() << std::endl;
        std::vector<int> len(probabilities.size());
	std::vector<double> P = probabilities;
	sort(P.begin(), P.end(), myfunction);
	std::vector<double> m_P = P;
        Huffman(C, len, P);
	P.clear();
	double coding_price=0;    
	for (int i = 0; i < C.size(); i++) {
		for (int z = 0; z < C.size(); z++) {
			if (probabilities[i] == m_P[z]) { // prob - исходный порядо букв в массиве,  m_P - отсортированный
				if (order[z] == -1) {
					order[z] = alphabet[i]; //создаем текущий порядок букв
					std::cout << alphabet[i] << " - ";
					for (int j = 0; j < C[z].size(); j++) {
						 std::cout << C[z][j];
					}
					 std::cout << "\n";


					coding_price += (m_P[z] * len[z]);
					P.push_back(m_P[z]);
					break;
				}
			}
		}
	}
	std::cout << "\n";
	std::cout << "Цена кодирования - " << coding_price << endl;
	std::cout << "Вектор последовательности символов - ";
	for (char val : order) {
            std::cout << val << ' ';
        }
        std::cout << '\n';   
	C_rectangular = transform_to_rectangle(C);
	int numRows = C_rectangular.size();
        int numCols = C_rectangular[0].size();
	for (const auto& row : C_rectangular) {
        for (int val : row) {
            std::cout << val << ' ';
        }
        std::cout << '\n';
    	}	
    	std::cout << std::endl;

	// for (int i = 1; i < world_size; ++i) {
	// MPI_Send(&numRows, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	// MPI_Send(&numCols, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	// for (int j = 0; j < numRows; ++j) {
	// 	MPI_Send(C_rectangular[j].data(), numCols, MPI_INT, i, 0, MPI_COMM_WORLD);
	// }
	// int orderSize = order.size();
	// MPI_Send(&orderSize, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	// MPI_Send(order.data(), orderSize, MPI_CHAR, i, 0, MPI_COMM_WORLD);
	// }
	// std::cout << "Отправилось с rank =0 "<< std::endl; 

	// CodingHuffman("Library.txt", "Coding", C_rectangular);

 //        // Принимаем закодированные строки от всех остальных узлов и записываем их в файл
 //        std::string encoded;
 //        for (int i = 1; i < world_size; ++i) {
 //            MPI_Status status;
 //            MPI_Probe(i, 0, MPI_COMM_WORLD, &status);

 //            int size;
 //            MPI_Get_count(&status, MPI_CHAR, &size);

 //            std::vector<char> received_data(size);
 //            MPI_Recv(received_data.data(), size, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

 //            std::string received_str(received_data.begin(), received_data.end());
 //            encoded += received_str;
        // }

        // std::ofstream output("Coding.txt");
        // output << encoded;
    }
    else {
        MPI_Send(symbols.data(), symbols.size(), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	std::cout << "Process " << world_rank << " sent message to process 0\n";
	    
	// // Принимаем информацию о таблице Хаффмана от процесса с rank = 0
	// MPI_Recv(&numRows, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	// MPI_Recv(&numCols, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	// C_rectangular.resize(numRows, std::vector<int>(numCols));
	// for (int i = 0; i < numRows; ++i) {
	//     MPI_Recv(C_rectangular[i].data(), numCols, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	// }
	
	// int orderSize;
	// MPI_Recv(&orderSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	// order.resize(orderSize);
	// MPI_Recv(order.data(), orderSize, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

 //        for (const auto& row : C_rectangular) {
 //        for (int val : row) {
 //            std::cout << val << ' ';
 //        }
 //        std::cout << '\n';
 //    	}	
 //    	std::cout << std::endl;	
     //std::string encoded = CodingHuffman("Library.txt", "Coding", C_rectangular);
     //MPI_Send(encoded.data(), encoded.size(), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }
	
//CodingRLE_MPI("Library.txt", "CodingRLE.txt");
//Coding RLE
{
    std::string file_RLE;
    if (world_rank == 0) {
        file_RLE = readFile("Library.txt");
    }

    int base_process, remainder;
    if (world_rank == 0) {
        int total_symbols = file_RLE.size();
        base_process = total_symbols / world_size;
        remainder = total_symbols % world_size;
    }

    MPI_Bcast(&base_process, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&remainder, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int start_symbol = world_rank * base_process + std::min(world_rank, remainder);
    int symbols_per_process = base_process + (world_rank < remainder ? 1 : 0);

    std::vector<char> chunk(symbols_per_process);

    if (world_rank == 0) {
        chunk.assign(file_RLE.begin() + start_symbol, file_RLE.begin() + start_symbol + symbols_per_process);
        for (int i = 1; i < world_size; i++) {
            int start_symbol_i = i * base_process + std::min(i, remainder);
            int symbols_per_process_i = base_process + (i < remainder ? 1 : 0);
            MPI_Status status;
            if (MPI_Send(file_RLE.data() + start_symbol_i, symbols_per_process_i, MPI_CHAR, i, 0, MPI_COMM_WORLD) != MPI_SUCCESS) {
                char err_string[MPI_MAX_ERROR_STRING];
                int len;
                MPI_Error_string(status.MPI_ERROR, err_string, &len);
                printf("MPI_Send failed with error: %s\n", err_string);
            }
        }
    } else {
        MPI_Status status;
        if (MPI_Recv(chunk.data(), chunk.size(), MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status) != MPI_SUCCESS) {
            char err_string[MPI_MAX_ERROR_STRING];
            int len;
            MPI_Error_string(status.MPI_ERROR, err_string, &len);
            printf("MPI_Recv failed with error: %s\n", err_string);
        }
    }

    std::cout << world_rank << ":   " << std::string(chunk.begin(), chunk.end()) << std::endl;
}
	
 if (world_rank == 0) {
	
 }
else {
	
}

    MPI_Finalize();
    return 0;
}
