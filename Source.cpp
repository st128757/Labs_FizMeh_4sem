#include <iostream>
#include <fstream>
#include <vector>
#include <string>

int main() {
	std::vector<double> tables;
	std::vector<double> S;

	std::string t1 = "C:\\programming\\Labs_FizMeh_4sem\\Stage-4.2\\1.txt";

	std::ifstream file(t1);
	if (!file.is_open()) {
		std::cerr << "Feilure " << t1 << "\n";
		std::cerr << "check path" << "\n";
		return 1;
	}
	std::cout << "Success. File is opened" << t1 << std::endl;

	double value;
	while (file >> value) {
		tables.push_back(value);
	
	}
	file.close();

	std::string t1S = "C:\\programming\\Labs_FizMeh_4sem\\Stage-4.2\\1S.txt";
	std::ofstream output_file(t1S);
	if (!output_file.is_open()) {
		std::cerr << "Failed to create output file: " << t1S << std::endl;
		return 1;
	}
	output_file << "sensibility of PVO\n";
	for (int i = 1; i < 6; ++i) {
		S.push_back(10 * i / (2 * sqrt(2) * tables[i - 1]));
		output_file << i << ") " << S[i - 1] << "\n";
	}
	output_file << "sensibility of PGO\n";
	for (int i = 1; i < 6; ++i) {
		S.push_back(10 * i / (2 * sqrt(2) * tables[i + 4]));
		output_file << i << ") " << S[i + 4] << "\n";
	}
	output_file << "max sensibility\n";
	for (int i = 1; i < 6; ++i) {
		S.push_back(10 * i / (2 * sqrt(2) * tables[i + 4]));
		output_file << i << ") " << S[i + 9] << "\n";
	}
	

	double S_x = 0.5 * (S[8] + S[9]);
	double S_y = 0.5 * (S[3] + S[2]);

	double U_x_mean = (18 + 22.5) * 0.5;
	double U_y_mean = (12.5 + 18.5) * 0.5;
	double U_x_err_rate = sqrt((pow(18 - U_x_mean, 2) + pow(22.5 - U_x_mean, 2)) / 2);
	double U_y_err_rate = sqrt((pow(12.5 - U_y_mean, 2) + pow(18.5 - U_y_mean, 2)) / 2);
	output_file << "Arithmetic means\n";
	output_file << "PGO: " << S_x << "\n" << "PVO: " << S_y << "\n";
	output_file << "Error rates of U\n";
	output_file << "PGO: " << U_x_err_rate << "\n" << "PVO: " << U_y_err_rate << "\n";

	double a[5] = { 0 };
	double K = a[0];
	for (int i = 0; i < 5; ++i) {
		a[i] = S[10 + i] / S[i];
		if (K < a[i]) {
			K = a[i];
		}
	}
	output_file << "koef max usilenya: " << K;
	output_file.close();
	return 0;
}