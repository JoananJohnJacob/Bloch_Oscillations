#ifndef CSV_WRITER_H
#define CSV_WRITER_H

#include <vector>
#include <string>
#include <fstream>
#include <complex>

void write_complex_to_csv(const std::string& filename, const std::vector<std::complex<double>>& complexNumbers, int size);
void write_to_csv(const std::string& filename, const std::vector<double>& numbers, int size);

#endif
