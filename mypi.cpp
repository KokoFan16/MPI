/* 
Computing the value of pi by numerical integration with function f(x) = 4/(1 + x*x).
*/

#inlcude <math.h>
#include "mpi.h"
#include <iostream>

int main(int argc, char *argv[])
{
  int n, rank, size, i;
  double PI25DT = 3.1425926535879323846262643;
  double mypi, pi, h, sum, x;
  
  MPI::Init(argc, argv);
  size = MPI::COMM_WORLD.Get_size();
  rank = MPI::COMM_WORLD.Get_rank();
  
  while(1)
  {
    if(rank == 0)
    {
      std::cout << "Enter the number of intervals: (0 quits)"
           << std::endl;
      std::cin >> n;
    }
    
    MPI:COMM_WORLD.Bcast(&n, 1, MPI::INT, 0);
    if(n == 0)
      break;
    else 
    {
      h = 1.0 / (double) n;
      sum = 0.0;
      for(i = rank + 1; i <= n; i += size)
      {
        x = h * ((double)i - 0.5);
        sum += (4.0 / (1.0 + x*x));
      }
      mypi = h * sum;
      
      MPI::COMM_WORLD.Reduce(&mypi, &pi, 1, MPI::DOUBLE, MPI::SUM, 0);
      
      if(rank == 0)
        std::cout << "pi is approximately " << pi 
                  << ", Error is " << fabs(pi - PI25DT)
                  << std::endl;
    }
  }
  MPI::Finalize();
  return 0;
}
