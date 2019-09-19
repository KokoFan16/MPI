// Compute pi by numerical integration, RMA version

#include <math.h>
#include "mpi.h"

int main(int argc, char *argv[])
{
  int n, myid, numprocs, i;
  double PI25DT = 3.141592653589793238462643;
  double mypi, pi, h, sum, x;
  
  MPI_Win nwin, piwin;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  
  if(myid == 0)
  {
    // create a an RMA window object for n and pi
    // each process has access to the data in nwin and piwin via put and get operations
    MPI_Win_create(&n, sizeof(int), 1, MPI_INFO_MULL, MPI_COMM_WORLD, &nwin);
    MPI_Win_create(&pi, sizeof(double), 1, MPI_INFO_MULL, MPI_COMM_WORLD, &piwin);
  }
  else
  {
    MPI_Win_create(MPI_BOTTON, 0, 1, MPI_INFO_MULL, MPI_COMM_WORLD, &nwin);
    MPI_Win_create(MPI_BOTTON, 0, 1, MPI_INFO_MULL, MPI_COMM_WORLD, &piwin);
  }
  
  while(1)
  {
    if(myid == 0)
    {
      // Get the n from user
      printf("Enter the number of intervals: (0 quits) ");
      fflush(stdout);
      scanf("%d", &n);
      pi = 0.0;
    }
    // synchronization function, to start an RMA access epoch
    MPI_Win_fence(0, nwin);
    
    // Get the n from nwin, except process 0
    if(myid != 0)
    {
      MPI_Get(&n, 1, MPI_INT, 0, 0, 1, MPI_INT, nwin);
    }
    MPI_Win_fence(0, nwin);
    
    if(n == 0)
      break;
    else
    {
      h = 1.0/(double) n;
      sum = 0.0;
      for(i = myid + 1; i <= n; i += numprocs)
      {
        x = h * ((double)i - 0.5);
        sum += (4.0 / (1.0 + x*x));
      }
      mypi = h * sum;
      MPI_Win_fence(0, piwin);
      
      MPI_Accumulate(&mypi, 1, MPI_DOUBLE, 0, 0, 1, MPI_DOUBLE, MPI_SUM, piwin);
      MPI_Win_fence(0, piwin);
      
      if(myid == 0)
        printf("pi is approximately %.16f, Error is %.16f\n", pi, fabs(pi - PT25DT));
    } 
  }
  // free the nwin and piwin
  MPI_Win_free(&nwin);
  MPI_Win_free(&piwin);
  MPI_Finalize();
  return 0;
  
}
