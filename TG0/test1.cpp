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
#include <cmath>

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

std::string DecodingRLE(const std::string& input_string, int& k) {
    k = input_string.size();
    std::string result;
    char symbol;
    int count;
    std::istringstream is(input_string);
    while (is >> symbol >> count) {
        for (int i = 0; i < count; i++) {
            if (symbol == 'q') result += '0';
            else if (symbol == 'w') result += '1';
            else if (symbol == 'e') result += '2';
            else if (symbol == 'r') result += '3';
            else if (symbol == 't') result += '4';
            else if (symbol == 'y') result += '5';
            else if (symbol == 'u') result += '6';
            else if (symbol == 'i') result += '7';
            else if (symbol == 'o') result += '8';
            else if (symbol == 'p') result += '9';
            else if (symbol == 'x') result += ' ';
            else result += symbol;
        }
    }
    return result;
}

std::string CodingHuffman(const std::string& s_output, const std::vector<std::vector<int>>& C, const std::string& substring) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::string result;
    for (char symbol : substring) {
        for (int i = 0; i < C.size(); ++i) {
            if (symbol == order[i]) {
                for (int j = 0; j < C[i].size(); ++j) {
                    if (C[i][j] != -1) {
                        result += std::to_string(C[i][j]);
                    }
                }
            }
        }
    }

    // Запись результата в файл
    std::ofstream output(s_output + "_part_" + std::to_string(rank)); // каждый процесс создает свой файл
    output << result;

    return result;
}

std::string DecodingHuffman(const std::string& s_input, const std::string& s_output, 
                            const std::vector<std::vector<int>>& C, int& k) {
    int rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Процессы открывают и читают соответствующий файл
    std::string file_content = readFile(s_input + "_part_" + std::to_string(rank));
    int total_symbols = file_content.size();

    std::string result;
    std::vector<int> code;

    for (int n = 0; n < total_symbols; ++n) {
        k += 1;  // Увеличиваем количество обработанных символов
        if (file_content[n] == '0') code.push_back(0);
        else if (file_content[n] == '1') code.push_back(1);

        for (int i = 0; i < C.size(); i++) {
            std::vector<int> temp = C[i];
            temp.erase(std::remove(temp.begin(), temp.end(), -1), temp.end()); // Удаляем все -1 из temp
            if (code == temp) {
                code.clear();
                result += order[i];  
                break;
            }
        }
    }
    std::ofstream output(s_output + "_part_" + std::to_string(rank)); // каждый процесс создает свой файл
    output << result;

    return result;
}

std::vector<char> readFile2(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    return std::vector<char>((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
}

bool Check(string s_input, string m_output) {
	ifstream input1(s_input);
	string coding; 
	ifstream input2(m_output);
	string decoding;
	bool flag = 1;
	while (getline(input1, coding) && getline(input2, decoding)) {
		int n = 0;
		while (n != coding.size()) {
			if (coding[n] != decoding[n]) flag = 0;
			n++;
		}
	}
	//system("pause");
	return flag;
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
    int numRows = 0;
    int numCols = 0;
    std::vector<std::vector<int>> C_rectangular;
    int k1=0, k2=0, k3 = 0, k4 = 0;
	
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

	for (int i = 1; i < world_size; ++i) {
	MPI_Send(&numRows, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	MPI_Send(&numCols, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	for (int j = 0; j < numRows; ++j) {
		MPI_Send(C_rectangular[j].data(), numCols, MPI_INT, i, 0, MPI_COMM_WORLD);
	}
	int orderSize = order.size();
	MPI_Send(&orderSize, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	MPI_Send(order.data(), orderSize, MPI_CHAR, i, 0, MPI_COMM_WORLD);
	}
	std::cout << "Отправилось с rank =0 "<< std::endl; 
	
    }
    else {
        MPI_Send(symbols.data(), symbols.size(), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	std::cout << "Process " << world_rank << " sent message to process 0\n";
	//MPI_Barrier(MPI_COMM_WORLD);
	    
	// Принимаем информацию о таблице Хаффмана от процесса с rank = 0
	MPI_Recv(&numRows, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(&numCols, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	C_rectangular.resize(numRows, std::vector<int>(numCols));
	for (int i = 0; i < numRows; ++i) {
	    MPI_Recv(C_rectangular[i].data(), numCols, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	
	int orderSize;
	MPI_Recv(&orderSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	order.resize(orderSize);
	MPI_Recv(order.data(), orderSize, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

     //    for (const auto& row : C_rectangular) {
     //    for (int val : row) {
     //        std::cout << val << ' ';
     //    }
     //    std::cout << '\n';
    	// }	
    	// std::cout << std::endl;	
    }

	
MPI_Barrier(MPI_COMM_WORLD);
std::vector<char> file_content;
    if (world_rank == 0) {
        file_content = readFile2("Library.txt");
    }

    int total_symbols = file_content.size();

    // Рассылаем total_symbols всем узлам
    MPI_Bcast(&total_symbols, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Теперь рассылаем содержимое файла всем узлам
    file_content.resize(total_symbols);
    MPI_Bcast(file_content.data(), total_symbols, MPI_CHAR, 0, MPI_COMM_WORLD);

    int base_process = total_symbols / world_size;
    int remainder = total_symbols % world_size;

    // Определение начала и конца обработки каждого процесса
    int start_symbol = world_rank * base_process + std::min(world_rank, remainder);
    int symbols_per_process = base_process + (world_rank < remainder ? 1 : 0);

    // Субстрока для этого процесса
    std::string substring(file_content.begin() + start_symbol, file_content.begin() + start_symbol + symbols_per_process);

    // Здесь ваш код для создания C_rectangular

    if (world_rank == 0) {
        std::string encoded = CodingHuffman("Coding", C_rectangular, substring);
        for (int i = 1; i < world_size; ++i) {
            MPI_Status status;
            MPI_Probe(i, 0, MPI_COMM_WORLD, &status);

            int size;
            MPI_Get_count(&status, MPI_CHAR, &size);

            std::vector<char> received_data(size);
            MPI_Recv(received_data.data(), size, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            std::string received_str(received_data.begin(), received_data.end());
            encoded += received_str;
        }

        std::ofstream output("Coding.txt");
        output << encoded;
    }
    else {
        std::string encoded = CodingHuffman("Coding", C_rectangular, substring);
        MPI_Send(encoded.data(), encoded.size(), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }
	
MPI_Barrier(MPI_COMM_WORLD);
 if (world_rank == 0) {
	std::string result1 = DecodingHuffman("Coding", "Decoding", C_rectangular, k1);
	std::vector<std::string> results(world_size);
	results[0] = result1;
        for (int i = 1; i < world_size; i++) {
            MPI_Status status;
            int result_size;
            MPI_Recv(&result_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
            char* buf = new char[result_size + 1];
            MPI_Recv(buf, result_size, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
            buf[result_size] = '\0';
            results[i] = std::string(buf);
            delete[] buf;
	}
	// Объединяем результаты и записываем их в один файл
        std::string combined_result_dec1;
        for (const auto& res : results) {
             combined_result_dec1 += res;
        }

        std::ofstream output_file("Decoding.txt");
        output_file <<  combined_result_dec1;
        output_file.close();
 }
else {
	std::string result1 = DecodingHuffman("Coding", "Decoding", C_rectangular, k1);
	int result_size = result1.size();
        MPI_Send(&result_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(result1.c_str(), result_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	std::cout << result1.size() ;
}

if (world_rank == 0) {
	bool flag1 = Check("Library.txt", "Decoding.txt");
	cout << "CheckHuffman:  " << flag1 << endl;
}

//////////////////////////////////////////////////////////////////////////////////////////// CodingRLE
if (world_rank == 0) {
	std::string encoded = CodingRLE(substring);
	std::ofstream output("CodingRLE.txt");
        output << encoded;
	std::ofstream output1("CodingRLE_part_" + std::to_string(world_rank) + ".txt");
        output1 << encoded;
	for (int i = 1; i < world_size; ++i) {
            MPI_Status status;
            MPI_Probe(i, 0, MPI_COMM_WORLD, &status);

            int size;
            MPI_Get_count(&status, MPI_CHAR, &size);

            std::vector<char> received_data(size);
            MPI_Recv(received_data.data(), size, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            std::string received_str(received_data.begin(), received_data.end());
            output << received_str;
        }
}
else{
	std::string encoded = CodingRLE(substring);
	MPI_Send(encoded.data(), encoded.size(), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	std::ofstream output1("CodingRLE_part_" + std::to_string(world_rank) + ".txt");
        output1 << encoded;
}

/////////////////////////////////////////////////////////////////////////////////////////////// DecodingRLE
    std::string encodedRLE;
    {
        std::ifstream input("CodingRLE_part_" + std::to_string(world_rank) + ".txt");
        std::stringstream ss;
        ss << input.rdbuf();
        encodedRLE = ss.str();
    }
    std::string decodedRLE = DecodingRLE(encodedRLE,k2);
    if (world_rank == 0) {
        std::ofstream output("DecodingRLE.txt");
        output << decodedRLE;
        for (int i = 1; i < world_size; ++i) {
            MPI_Status status;
            MPI_Probe(i, 0, MPI_COMM_WORLD, &status);

            int size;
            MPI_Get_count(&status, MPI_CHAR, &size);

            std::vector<char> received_data(size);
            MPI_Recv(received_data.data(), size, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            std::string received_str(received_data.begin(), received_data.end());
            output << received_str;
        }
    } else {
        MPI_Send(decodedRLE.data(), decodedRLE.size(), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }
	bool flag2 = Check("Library.txt", "DecodingRLE.txt" );
	std::cout << "CheckRLE:  " << flag2 << std::endl;
	std::cout << "Сompression ratio RLE: " << 10000. / k2 << std::endl;
	std::cout << "Сompression ratio : " << 10000. / k1 << std::endl;
    MPI_Finalize();
    return 0;
}

















