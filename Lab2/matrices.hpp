

#ifndef matrices_hpp
#define matrices_hpp

#include <stdio.h>
#include <vector>
#include <iomanip>
#include <iostream>


std::vector<std::vector<float>> RandomMatrix(unsigned n);
std::vector<std::vector<float>> HilbertMatrix(unsigned n);
std::vector<std::vector<float>> RandomMatrix2(unsigned n);
std::vector<std::vector<float>> HilbertMatrix2(unsigned n);

void printMatrix(std::vector<std::vector<float>> matrix, unsigned n);

#endif /* matrices_hpp */
