#include "matrices.hpp"

std::vector<std::vector<float>> RandomMatrix(unsigned n)
{
    std::vector<std::vector<float>> a (n, std::vector<float>(n, 0));
    srand( static_cast<unsigned int>(time(nullptr)));
    
    for( int i = 0; i < n; ++i)
        {for( int j = 0;  j < n; ++j)
           {a[i][j] = rand()%10;}
        }
       
    return a;
    
}


std::vector<std::vector<float>> RandomMatrix2(unsigned n) //returns random matrix with dominant diagonal
{
    std::vector<std::vector<float>> a (n, std::vector<float>(n, 0));
    srand( static_cast<unsigned int>(time(nullptr)));
    
    for( int i = 0; i < n; ++i)
    {
        {for( int j = 0;  j < n; ++j)
           {a[i][j] = rand()%10;}
        }
    }
    for( int i = 0; i < n; ++i)
    {
        float sum = 0;
        for( int j = 0;  j < n; ++j)
           {
               if(i!=j)
               {
                   sum+=a[i][j];
               }
           }
        if (sum>a[i][i])
        {
            a[i][i] = sum;
        }
        
    }
    
    
    return a;
    
}

std::vector<std::vector<float>> HilbertMatrix(unsigned n)
{
  std::vector<std::vector<float>> a (n, std::vector<float>(n, 0));
  for(int i = 0; i<n; i++)
    {
        for(int j = 0; j <n; j++)
        {
            a[i][j]= (float)1.0 /
                ((i + 1) + (j + 1) - 1.0);
        }
    }
    
  return a;
   
}
 
std::vector<std::vector<float>> HilbertMatrix2(unsigned n)//returns hilbert matrix with dominant diagonal
{
  std::vector<std::vector<float>> a (n, std::vector<float>(n, 0));
  for(int i = 0; i<n; i++)
    {
        for(int j = 0; j <n; j++)
        {
            a[i][j]= (float)1.0 /
                ((i + 1) + (j + 1) - 1.0);
        }
    }
    for( int i = 0; i < n; ++i)
    {
        float sum = 0;
        for( int j = 0;  j < n; ++j)
           {
               if(i!=j)
               {
                   sum+=a[i][j];
               }
           }
        if (sum>a[i][i])
        {
            a[i][i] = sum;
        }
        
    }
  return a;
   
}
 

void printMatrix(std::vector<std::vector<float>> matrix, unsigned n)
{
    for(int i = 0; i<n; i++)
    {
        for(int j = 0; j <n; j++)
        {
            std::cout <<std::setfill(' ')<<std::setw(12);
            std::cout<<std::setprecision(6)<<matrix[i][j];
           
        }
        std::cout<<"\n";
    }
}

