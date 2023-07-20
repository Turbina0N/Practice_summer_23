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
static char order[19] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
static char orderNew[19] = { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 };

vector<double> CreateFile() {
	FILE* pFile;
	fopen_s(&pFile, "Library.txt", "w");
	int i;
	vector<int> entry(19, 0);
	vector<double> probability(19, 0);
	char num;
	int index;
	for (i = 0; i < N; ++i)
	{
		index = rand() % 19 ;
		num = arr[index];
		entry[index]++;
		fwrite(&num, 1, sizeof(num), pFile);
	}
	for (int i = 0; i < probability.size(); i++) {
		probability[i] = entry[i] / (double)N;
	}

	fclose(pFile);
	return probability;
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
	return j; //номер вставляемого элемента
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
	string string; 
	ofstream output(s_output);
	
	while (getline(input, string)) { 
		int n = 0;
		while (n != string.size()) {
			for (int i = 0; i < C.size(); i++) {
				if (string[n] == order[i]) {
					for (int j = 0; j < C[i].size(); j++)
						result += to_string(C[i][j]);  
					n++;
				}
				
			}
		}
	}
	output << result;
	return result;
	//system("pause");
}

string CodingHuffmanNew(string s_input, string s_output, vector<vector<int>> C) {
	string result;
	ifstream input(s_input);
	string string;
	ofstream output(s_output);
	for (int i = 0; i < 19; i++) {
		if (order[i] == '0') orderNew[i] = 'q';
		else if (order[i] == '1') orderNew[i] = 'w';
		else if (order[i] == '2') orderNew[i] = 'e';
		else if (order[i] == '3') orderNew[i] = 'r';
		else if (order[i] == '4') orderNew[i] = 't';
		else if (order[i] == '5') orderNew[i] = 'y';
		else if (order[i] == '6') orderNew[i] = 'u';
		else if (order[i] == '7') orderNew[i] = 'i';
		else if (order[i] == '8') orderNew[i] = 'o';
		else if (order[i] == '9') orderNew[i] = 'p';
		else if (order[i] == ' ') orderNew[i] = 'x';
		else orderNew[i] = order[i];
	}
	while (getline(input, string)) {
		int n = 0;
		while (n != string.size()) {
			for (int i = 0; i < C.size(); i++) {
				if (string[n] == orderNew[i]||string[n] == order[i]) {
					for (int j = 0; j < C[i].size(); j++)
						result += to_string(C[i][j]);
					n++;
				}

			}
		}
	}
	output << result;
	return result;
	//system("pause");
}



//bool CheckHuffman() {
//	ifstream input1();
//	string coding; 
//	ifstream input2(); 
//	string decoding; 
//	bool flag = 1;
//	while (getline(input1, coding) && getline(input2, decoding)) {
//		int n = 0;
//		while (n != coding.size()) {
//			if (coding[n] != decoding[n]) flag = 0;
//			n++;
//		}
//	}
//	//system("pause");
//	return flag;
//}

string DecodingHuffman(string s_input, string s_output, vector<vector<int>> C, int& k) {
	string result;
	ifstream input(s_input);
	string string;
	ofstream output(s_output);
	vector<int> code;
	while (getline(input, string)) {
		k += string.size();
		int n = 0;
		bool flag = 0;
		int j = 0;
		while (n != string.size()) {
			if (string[n] == '0') code.push_back(0);
			if (string[n] == '1') code.push_back(1);
			for (int i = 0; i < C.size(); i++) {
				if (code == C[i]) {
					flag = 1;//нашелся
					j = i;
				}
			}
			if (flag) {
				code.clear();
				result += order[j];  //И записываешь эту строку в выходной файл
				flag = 0;
			}
			n++;
		}
	}
	output << result;
	return result;
	//system("pause");
}


string CodingRLE(string s_input, string s_output) {

	ifstream input(s_input);
	string str; 
	string allstr; 
	ofstream output(s_output);
	while (getline(input, str)) { 
		allstr += str;
	}

	string result;
	int count = 1;
	for (int i = 1; i < allstr.size(); i++) {
		if (allstr[i] == allstr[i - 1]) {
			count++;
		}
		else {
			if (allstr[i - 1] == '0') result += 'q';
			else if (allstr[i - 1] == '1') result += 'w';
			else if (allstr[i - 1] == '2') result += 'e';
			else if (allstr[i - 1] == '3') result += 'r';
			else if (allstr[i - 1] == '4') result += 't';
			else if (allstr[i - 1] == '5') result += 'y';
			else if (allstr[i - 1] == '6') result += 'u';
			else if (allstr[i - 1] == '7') result += 'i';
			else if (allstr[i - 1] == '8') result += 'o';
			else if (allstr[i - 1] == '9') result += 'p';
			else if (allstr[i - 1] == ' ') result += 'x';
			else result += allstr[i - 1];
			result += to_string(count);
			count = 1;
		}
	}
	if (allstr.back() == '0') result += 'q';
	else if (allstr.back() == '1') result += 'w';
	else if (allstr.back() == '2') result += 'e';
	else if (allstr.back() == '3') result += 'r';
	else if (allstr.back() == '4') result += 't';
	else if (allstr.back() == '5') result += 'y';
	else if (allstr.back() == '6') result += 'u';
	else if (allstr.back() == '7') result += 'i';
	else if (allstr.back() == '8') result += 'o';
	else if (allstr.back() == '9') result += 'p';
	else if (allstr.back() == ' ') result += 'x';
	else result += allstr.back();
	result += to_string(count);


	output << result;
	return result;
}

string DecodingRLE(string s_input, string m_output, int& k) {
	ifstream input(s_input);
	string str; 
	string allstr; 
	ofstream output(m_output);

	while (getline(input, str)) { 
		allstr += str;
	}
	k = allstr.size();

	string result;
	char symbol;
	int count;
    istringstream is(allstr);
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

	output << result;
	return result;

	/*char symbol;
	int count;
	string num;
	int i = 0;
	while (i != allstr.size()) {
	symbol = allstr[i];
		i++;
		while (allstr[i] != '*') {
		num += allstr[i];
		i++;
	}
	count = stoi(num);
	for (int t = 0; t < count; t++) {
		output << symbol;
	}
	num.clear();
	i++;
}*/

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

string Twice1(int& k) {
	//string strHuffman = CodingHuffman("Library.txt", "CodingTwice1.txt", C);
	string strRLE = CodingRLE("Coding.txt", "ResultCodingTwice1");
	ifstream input("ResultCodingTwice1");
	string str;
	int num = 0;
	while (getline(input, str)) {
		num += str.size();
	}
	k = num;
	return strRLE;
}

string Twice2(int& k, vector<vector<int>> C) {
	string strHuffman = CodingHuffmanNew("CodingRLE.txt", "ResultCodingTwice2.txt", C);
	//string strRLE = CodingRLE("Coding.txt", "ResultCodingTwice1");
	ifstream input("ResultCodingTwice2.txt");
	string str;
	int num = 0;
	while (getline(input, str)) {
		num += str.size();
	}
	k = num;
	return strHuffman;
}

int main() {
	SetConsoleCP(1251);         // !!!
	SetConsoleOutputCP(1251);   // !!!
	//setlocale(LC_ALL, "Russian");
	srand(time(NULL));
	vector<double> probability = CreateFile();
	double coding_price=0; // цена кодирования
	vector<vector<int>> C;
	vector<int> len;
	//vector<double> P = {0.2, 0.2, 0.19, 0.12, 0.11, 0.09, 0.09}; //проверка алгоритма Хаффмена
	vector<double> P = probability;
	//j = UP(P, 0.1);

	sort(P.begin(), P.end(), myfunction);
	vector<double> m_P = P;

	C.resize(P.size());
	len.resize(P.size());
	Huffman(C, len, P);

	/*for (int i = 0; i < C.size(); i++) {
				for (int j = 0; j < C[i].size(); j++) {
					cout << C[i][j];
				}
				cout << "\n";
	}*/
	P.clear();

	for (int i = 0; i < C.size(); i++) {
		for (int z = 0; z < C.size(); z++) {
			if (probability[i] == m_P[z]) { // prob - исходный порядо букв в массиве,  m_P - отсортированный
				if (order[z] == -1) {
					order[z] = arr[i]; //создаем текущий порядок букв
					cout << arr[i] << " - ";
					for (int j = 0; j < C[z].size(); j++) {
						cout << C[z][j];
					}
					cout << "\n";


					coding_price += (m_P[z] * len[z]);
					P.push_back(m_P[z]);
					break;
				}
			}
		}
	}

	cout << "\n";
	cout << "Цена кодирования - " << coding_price << endl;

	//Новый порядок с их кодированием
	//for (int i = 0; i < C.size(); i++) {
	//	cout << order[i] << " - ";
	//	for (int j = 0; j < C[i].size(); j++) {
	//		cout << C[i][j];
	//	}
	//	cout << "\n";
	//}

	int k1=0, k2=0, k3 = 0, k4 = 0;

	string codingHuffman = CodingHuffman("Library.txt", "Coding.txt", C);
	string decodingHuffman = DecodingHuffman("Coding.txt", "Decoding.txt", C, k1);
	bool flag1 = Check("Library.txt", "Decoding.txt");
	cout << "CheckHuffman:  " << flag1 << endl;

	string codingRLE = CodingRLE("Library.txt", "CodingRLE.txt");
	string decodingRLE = DecodingRLE("CodingRLE.txt", "DecodingRLE.txt", k2);
	bool flag2 = Check("Library.txt", "DecodingRLE.txt" );
	cout << "CheckRLE:  " << flag2 << endl;
	cout << "Сompression ratio RLE: " << 10000. / k2 << endl;
	cout << "Сompression ratio : " << 10000. / k1 << endl;

	string v1 = Twice1(k3);
	string v2 = Twice2(k4, C);

	cout << "Сompression ratio : " << 10000. / k3 << endl;
	cout << "Сompression ratio : " << 10000. / k4 << endl;

	return 0;
}