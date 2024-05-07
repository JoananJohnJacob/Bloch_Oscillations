#include "csv_writer.h"
#include <iostream>
#include <fstream>

void write_complex_to_csv(const std::string& filename, const std::vector<std::complex<double>>& complexNumbers, int size) {
    std::ofstream outputFile(filename);
    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to Open file" << filename;
        return;
    }
    for (const auto& num : complexNumbers) {
        outputFile << num.real() << "," << num.imag() << std::endl;
    }
    outputFile.close();
}

void write_to_csv(const std::string& filename, const std::vector<double>& numbers, int size) {
    std::ofstream outputFile(filename);
    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to Open file" << filename;
        return;
    }
    for (const auto& num : numbers) {
        outputFile << num << std::endl;
    }
    outputFile.close();
}
