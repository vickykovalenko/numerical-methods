#include "functions.hpp"
#include "matrices.hpp"
#define eps 0.001

void GaussElim() //Gaussian Elimination wuth LU decomposition and partial pivoting
{   int n;
    std::cout<<"Enter quantity of equations:\t";
    std::cin>>n;
    int l,m,p;
    
    std::vector<std::vector<float>> a (0, std::vector<float>(0));
    std::vector<std::vector<float>> c (0, std::vector<float>(0));
    std::vector<std::vector<float>> U(0, std::vector<float>(0));
    std::vector<std::vector<float>> L(0, std::vector<float>(0));
    
    std::vector<float> b, d, q;
    std::vector<float> x(n, 0);
    float sum;
    /*for (int i=0; i<n ; i++)
    {
       std::vector<float> row;
       std::cout <<"\nEnter Row"<<i<<" coefficients:\n";
        for(int j = 0; j<n; j++)
        {
            
            std::cout<<i<<j<<"=";
             float num;
            std::cin>>num;
            row.push_back(num);
            
        }
        a.push_back(row);
        c=a;
        std::cout<<"_______________________________________________________________";
    }*/
    std::cout<<"Choose type of matrix: 1 - Random, 2 - Hilbert"<<std::endl;
    int choice;
    std::cin>>choice;
    if(choice == 1)
    {
         std::vector<std::vector<float>> matrix = RandomMatrix(n);
         printMatrix(matrix, n);
         a = matrix;
         c=a;

    }
    else if(choice == 2)
    {
        std::vector<std::vector<float>> matrix = HilbertMatrix(n);
        printMatrix(matrix, n);
        a = matrix;
        c=a;

    }
    else{
        std::cout<<"OOps!Try again!"<<std::endl;
    }

    std::cout<<"\nEnter b (the right side of the equations):\n";
    for (int i=0 ; i<n ; i++)
    {
        std::cout<<"\tb"<<i<<"=";
        float num;
        std::cin>>num;
        b.push_back(num);
        q=b;
    }
    std::cout<<"_______________________________________________________________";
    
    //pivoting
    std::cout<<"\npivoting"<<std::endl;
    
    for (int i=0 ; i<n ; i++){
        m = abs(a[i][i]);
        for (int j=i ; j<n ; j++){
            l=abs(a[j][i]);
            if(l>m)
            {
                std::cout<<"m = "<<m<<std::endl;
                std::cout<<"l = " <<l<<std::endl;
                m=l;
                for (int k=0 ; k<n ; k++){
                    a[i][k]=c[j][k];
                    a[j][k]=c[i][k];
                    c[i][k]=a[i][k];
                    c[j][k]=a[j][k];
                }
                
                b[i]=q[j];
                b[j]=q[i];
                q[i]=b[i];
                q[j]=b[j];
            }
        }
    }
    std::cout<<"\n\nThe matrix a after pivoting is:";
    for (int i=0 ; i<n ; i++){
        std::cout<<"\n";
        for (int j=0 ; j<n ; j++)
         {
             std::cout<<"\ta"<<i<<j<<"="<<a[i][j];
         }
    }
    
     std::cout<<"\nThe b matrix after pivoting is:";
     for (int i=0 ; i<n ; i++)
     { std::cout<<"\n\tb"<<i<<"="<<b[i];}
     std::cout<<"\n_______________________________________________________________";
    
    //"LU decomposition"
    //L & U Initial Values //U[n][n]:
    for (int i=0 ; i<n ; i++)
    {
        std::vector<float> row;
        for (int j=0; j<n ; j++)
        {
           
            if (i==0){
               
                float num = a[i][j];
                
                //U[i][j]=a[i][j];
                row.push_back(num);
                
            }
            else
            {
                 row.push_back(0);
                //U[i][j]=0;
            }
          
        }
        U.push_back(row);
    }
    //L[n][n]:
    for (int i=0; i<n ; i++)
    {
        std::vector<float> row2;

        for (int j=0; j<n ; j++){
        if(i==j)
        {//L[i][j]=1;
            row2.push_back(1);
        }
        else{
            //L[i][j]=0;
            row2.push_back(0);
        }
        }
         L.push_back(row2);
    }
    
    //" LU-decomposition"
    for (int i=1 ; i<n ; i++){
        for (int j=0 ; j<n ; j++){
            sum = 0;
            if(i>j) {
                for (int k=0 ; k<n ; k++){   //arrays in row i
                    sum = sum+L[i][k]*U[k][j];
                }
                L[i][j]=(a[i][j]-sum)/U[j][j];
                }
                else {
                    sum=0;
                    for(int k = 0; k<n; k++){
                        sum = sum+L[i][k]*U[k][j];
                    }
                    U[i][j]=a[i][j]-sum;
                }
        }
    }
    std::cout<<"\n\nL Matrix is:";
    for (int i=0 ; i<n ; i++){
       std::cout<<"\n";
        for (int j=0 ; j<n ; j++){
            std::cout<<"\tL"<<i<<j<<"="<<L[i][j];
            }
        }
    std::cout<<"\n\nU Matrix is:";
    for (int i=0 ; i<n ; i++){
        std::cout<<"\n";
        for (int j=0 ; j<n ; j++){
            std::cout<<"\tU"<<i<<j<<"="<<U[i][j];
            }
        }
    std::cout<<"_______________________________________________________________";
    p=0;
    for (int i=1 ; i<n ; i++)
    {
        if(U[i][i]==0)
        {
        p=p+1; }
    }
    if(p>0) {
    std::cout<<"\nSystem has infinitely many solutions";
    }
    else{
        std::cout<<"\n_______________________________________________________________";
        //calculating d for each line(Ld=b)

        
        float num = b[0];
        d.push_back(num);
        
        for (int i = 1 ; i<n ; i++){
            sum=0;
            for (int j=0 ; j<=i-1 ; j++){
             sum = sum+d[j]*L[i][j];
            }
             float num=b[i]-sum;
            d.push_back(num);
        }
        std::cout<<"\n\nThe d values are:";
        for (int i=0 ; i<n ; i++){
            std::cout<<"\n\td"<<i<<"="<<d[i];
       
        }
        std::cout<<"\n_______________________________________________________________";
        //calculating x(Ux=d)
        //num = d[n]/U[n][n];
        x[n-1]=d[n-1]/U[n-1][n-1];
        for (int i=n-2; i>=0 ; i--)
        {
            sum=0;
        for(int j=i+1 ; j<n ; j++){
            sum=sum+U[i][j]*x[j];
            }
            x[i]=(d[i]-sum)/U[i][i];
        }
        std::cout<<"\n\nThe x values are:";
        for (int i=0 ; i<n ; i++)
        {  std::cout<<"x"<<i<<"="<<x[i]<<std::endl;
        }
    }
}
void Jacobi()
{
    unsigned n;
    
    std::cout << "\nEnter quantity of equations:\nn = ";
    std::cin >> n;
    
    std::vector<float> x(n);
    std::vector<std::vector<float>> a (0, std::vector<float>(0));
    std::vector<float> b(n);
   
    std::cout<<"Choose type of matrix: 1 - Random, 2 - Hilbert"<<std::endl;
    int choice;
    std::cin>>choice;
    
    if(choice == 1)
    {
        std::vector<std::vector<float>> matrix = RandomMatrix2(n);
        printMatrix(matrix, n);
        a = matrix;
    }
    else if(choice == 2)
    {
        std::vector<std::vector<float>> matrix = HilbertMatrix2(n);
        printMatrix(matrix, n);
        a = matrix;
    }
    else{
           std::cout<<"OOps!Try again!"<<std::endl;
    }
    std::cout << "\nEnter b-vector\n";
    for (auto& b_elem : b)
    {
        std::cin >> b_elem;
    }
    float allowed_error = 0.01;
    
    //solve(a, b, x, allowed_error);
    const unsigned n2 = x.size();
    std::vector<float> tmp_x(n);
    
    float error;
    
    do
    {
        error = 0;
        
        tmp_x = b;
        for (unsigned i = 0; i < n2; ++i)
        {
            for (unsigned j = 0; j < n2; ++j)
            {
                if (i != j)
                {
                    tmp_x[i] -= a[i][j] * x[j];
                }
            }
            
            //оновити x[i] та порахувати похибку
            const float x_updated = tmp_x[i] / a[i][i];
            const float e = fabs(x[i] - x_updated);
            x[i] = x_updated;
            if (e > error) { error = e; }
        }
    }
    while (error > allowed_error);
    
    std::cout << "\nРозв'язок системи:: \n";
    for (unsigned i = 0; i < n; ++i)
    {
        std::cout << "x[" << i << "]=  " << x[i]<<std::endl;
    }
}
bool converge(std::vector<float> xk, std::vector<float> xkp, int n)
{
    bool b = true;
    for (int i = 0; i < n; i++)
    {
        if (abs(xk[i]-xkp[i]) > eps)
        {
            b = false;
            break;
        }
    }
    return b;
}
void GaussSeidel()
{
    std::vector<std::vector<float>> a (0, std::vector<float>(0));
    std::vector<float> b;

    int n = 0;
    std::cout << "Enter number of equations : ";
    std::cin >> n;
    std::vector<float> x(n, 0);
    std::vector<float> m(n, 0);

    std::cout<<"Choose type of matrix: 1 - Random, 2 - Hilbert"<<std::endl;
    int choice;
    std::cin>>choice;
    
    if(choice == 1)
    {
        std::vector<std::vector<float>> matrix = RandomMatrix2(n);
        printMatrix(matrix, n);
        a = matrix;
    }
    else if(choice == 2)
    {
        std::vector<std::vector<float>> matrix = HilbertMatrix2(n);
        printMatrix(matrix, n);
        a = matrix;
    }
    else{
           std::cout<<"OOps!Try again!"<<std::endl;
    }
    std::cout<<"\nEnter b (the right side of the equations):\n";
    for (int i=0 ; i<n ; i++)
    {
        std::cout<<"\tb"<<i<<"=";
        float num;
        std::cin>>num;
        b.push_back(num);
    }
    float y;
    int q, p ,i, j;
        for (int i=0;i<n;i++)                //Loop that calculates x1,x2,...xn
        {
           std::cout << "\nEnter the no. of iteration : ";
             std::cin >> q;
             while (q> 0) {
                for (i = 0; i < n; i++) {
                   x[i] = (b[i] / a[i][i]);
                   for (int j = 0; j < n; j++) {
                      if (j == i)
                         continue;
                      x[i] = x[i] - ((a[i][j] / a[i][i]) * m[j]);
                      m[i] = x[i];
                   }
                   std::cout<<"x"<<i + 1 << "="<<x[i]<<" ";
                }
                std::cout << "\n";
                q--;
             }
        }
    }

