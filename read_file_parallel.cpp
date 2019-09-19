#include "mpi.h"
#include <stdio.h>

int main(int argc, char *argv[])
{
  int myrank, numprocs, bufsize, *buf, count;
  NPI_File thefile;
  MPI_Status status;
  MPI_Offset filesize;
  
  MPI_Init(&argc, &argv); // Init the mpi environment
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank); // Get the current rank of process
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs); // Get the number of processes
  MPI_File_Open(MPI_COMM_WORLD, "testfile", MPI_MODE_RDONLY, MPI_INFO_NULL, &thefile); //open a file
  MPI_File_get_size(thefile, &filesize); // get the size of file in bytes
  filesize = filesize/sizeof(int); //in number of int
  bufsize = filesize/numprocs + 1; //local number to read
  buf = (int *) malloc(bufsize * sizeof(int)); // local buffer
  MPI_File_set_view(thefile, myrank * bufsize * sizeof(int),  //set the individual read pointer for each process
                   MPI_INT, MPI_INT, "native", MPI_INFO_NULL);
  MPI_File_read(thefile, buf, bufsize, MPI_INT, &status); // read file in parallel
  MPI_Get_count(&status, MPI_INT, &count); // get the the number of read ints
  printf("process %d read %d ints\n", myrank, count);
  MPI_File_close(&thefile); //close file
  MPI_Finalize(); //end mpi
  return 0;  
}
