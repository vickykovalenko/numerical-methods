#include <iostream>
#include "matrices.hpp"
#include "functions.hpp"

int main() {
    int choice;
    std::cout<<"Choose an algorithm: 1 - Gaussian elimination, 2 - Jacobi method, 3 - Gauss-Seidel method"<<std::endl;
    std::cin>>choice;
    if (choice == 1)
    {GaussElim();}
    else if(choice == 2)
    {Jacobi();}
    else if(choice == 3)
    {GaussSeidel();}
    else {
        std::cout<<"OOps!Try again!"<<std::endl;
    }
}
