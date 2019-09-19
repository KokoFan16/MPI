#include "mpi.h"

#define FILESIZE (1024 * 1024)

int main(int argc, char **argv)
{
  int *buf, rank, nprocs, nints, bufsize;
  MPI_File fh;
  MPI_Status status;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  
  bufsize = FILESIZE/nprocs; // int bytes
  buf = (int *) malloc(bufsize);
  nints = bufsize/sizeof(int); // the number of ints
  offset = rank * bufsize; // offset for each process
  
  MPI_File_open(MPI_COMM_WORLD, "testfile", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  // offset is passed as an argument to MPI_File_read_at
  MPI_File_read_at(fh, offset, buf, nints, MPI_INT, &status);
  MPI_File_close(&fh);
  
  free(buf);
  MPI_Finalize();
  return 0;
