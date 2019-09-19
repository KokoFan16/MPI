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
  
  MPI_File_open(MPI_COMM_WORLD, "testfile", MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
  // Set the position of individual pointer
  MPI_File_seek(fh, rank * bufsize, MPI_SEEK_SET);
  MPI_File_write(fh, buf, nints, MPI_INT, &status);
  MPI_File_close(&fh);
  
  free(buf);
  MPI_Finalize();
  return 0;
