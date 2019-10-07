#include <unistd.h>
#include <stdarg.h>
#include <stdint.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <time.h>
#include <string.h>
#include <mpi.h>

// Parse arguments
static void parse_args(int argc, char **argv);
// Global dimentsion of picture
static int gbs[3];
// Name of file to read
static char file_name[512];
// Name of file to write
static char write_file_name[512];
// Image format e.g., uint8
static int image_format;
// The matrix size which is used to filter the image
static int matrix_size;
// The size of each pixel. e.g., uint8 is 1 byte.
static int dt_size;
// The border of the image based on the size of filtering matrix. e.g., matrix 3*3*3, the border is 1
static int border;





// Generation a filtering matrix based on the parameter -m (matrix_size)
// j, i, k are cloumn, row and slice respectively.
// Return a array[matrix_size * matrix_size * matrix_size].
// The values in the array are the data used in filtering which is around the input data.
unsigned char * generate_filtering_matrix(int k, int i, int j, int dt_size, unsigned char *read_buf)
{
    
    typedef union {
        char c;
        short s;
        int i;
        float f;
    } datatype;
    
    
    datatype data1;
    datatype data2;
    
    
    int mid = matrix_size/2;
    int avg_value = matrix_size * matrix_size * matrix_size;
    
    // Return value
    //int *ret = malloc(avg_value * sizeof(int));
    unsigned char * ret = malloc(dt_size);

    if(dt_size == 2)
    {
        data1. s = 0;
        data2.s = 0;
    }
    
    int index = 0;

    int a = 0, b = 0, c = 0;
    for(a = mid*(-1); a < mid + 1; a++)
    {
        for(b = mid*(-1); b < mid + 1; b++)
        {
            for(c = mid*(-1); c < mid + 1; c++)
            {
                //buf[(k + 1) * gbs[0] * gbs[1] + i * gbs[1] + j]
                //ret[index] = (k + a) * gbs[0] * gbs[1] + ((i + b) * gbs[1]) + (j + c);
                //index += 1;
                
                index = (k + a) * gbs[0] * gbs[1] + ((i + b) * gbs[1]) + (j + c);
                
                if(dt_size == 1)
                {
                    ret[0] += read_buf[index]/avg_value;
                }
                
                if(dt_size == 2)
                {
                    memcpy(&data1.s, &read_buf[index * dt_size], dt_size);
                    data2.s += data1.s/avg_value;
                }
            }
        }
    }
    
    if(dt_size == 2)
    {
        memcpy(ret, &data2.s, dt_size);
    }
    
    return ret;
}


// Main funciton
int main(int argc, char **argv)
{
    
    parse_args(argc, argv);
    
    
    // This is the number of bytes per pixel, e.g., uint16 is 2
    dt_size = image_format/8;
    border = matrix_size/2;
    
    // Initialize the MPI execution environment
    MPI_Init(&argc, &argv);
    
    // rank is rank of the calling process in the group of comm (integer)
    // nprocs is number of processes in the group of comm (integer)
    int rank, nprocs, bufsize, n, ncells;
    int avg_value = matrix_size * matrix_size * matrix_size;
    
    // Local buffer for each process
    unsigned char * buf;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    // READ THE FILE IN PARALLEL
    MPI_File fh[2];
    MPI_Status status[2];
    
    // Returns an elapsed time on the calling processor
    double io_read_start = MPI_Wtime();
    
    // The local buffer size for each process.
    bufsize = (gbs[0] * gbs[1]) * (gbs[2]/nprocs);
    // The remaining slices which cannot be divisble
    int overflow = gbs[2]%nprocs;
    
    // The last process is used to deal with the overflow part
    if(rank == nprocs - 1)
    {
        //bufsize + overflow * gbs[0]*gbs[1]) * dt_size
        buf = malloc((bufsize + overflow * gbs[0] * gbs[1] + 2 * border * gbs[0] * gbs[1]) * dt_size);
        n = (bufsize + overflow*gbs[0]*gbs[1]) * dt_size;
    }
    else
    {
        // The extra space (gbs[0] * gbs[1] * 2 * border) is allocated for ghost cells
        // buf[0] ~ buf[gbs[0] * gbs[1] * border)] is for the ghost cells from last rank
        // buf[bufsize + gbs[1] * gbs[2] * border)] ~ buf[bufsize + 2 * gbs[1] * gbs[2] * border)] is for the ghost cells from the next rank
        buf = malloc((bufsize + 2 * border * gbs[0]*gbs[1]) * dt_size);
        // The number of short to be read
        n = bufsize * dt_size;
    }
    
    //Open file in parallel
    MPI_File_open(MPI_COMM_WORLD, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh[0]);
    // Set the read location for each rank
    MPI_File_seek(fh[0], rank * bufsize * dt_size, MPI_SEEK_SET);
    // Store the data from buf[border * gbs[0] * gbs[1]], since the first border slice leaves for ghost cells
    MPI_File_read(fh[0], &buf[border * gbs[0] * gbs[1] * dt_size], n, MPI_UNSIGNED_CHAR, &status[0]);
    
    MPI_File_close(&fh[0]);
    
    
    double io_read_end = MPI_Wtime();
    
    
    // GHOST (HALO) CELL EXCHANGE
    double communication_start = MPI_Wtime();
    
    MPI_Request reqs[4];
    MPI_Status stats[4];
    
    // The number of ghost cells for each ghost block
    ncells = border * gbs[0] * gbs[1] * dt_size;
   
    // MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
    // MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request)
    // Start ghost exchange
    if (rank == 0)
    {
        // For the first block, it just exchanges the last broder slices, e.g., for uint16, exchange last 2 slices
        MPI_Isend(&buf[bufsize * dt_size], ncells, MPI_UNSIGNED_CHAR, 1, 0, MPI_COMM_WORLD, &reqs[0]);
        MPI_Irecv(&buf[(bufsize + border * gbs[0]*gbs[1]) * dt_size], ncells, MPI_UNSIGNED_CHAR, 1, 1, MPI_COMM_WORLD, &reqs[1]);
        
        MPI_Wait(&reqs[0], &stats[0]);
        MPI_Wait(&reqs[1], &stats[1]);
    }
    
    else if (rank == (nprocs - 1))
    {
        // For the last block, it just exchanges the first broder slices
        MPI_Isend(&buf[border * gbs[0]*gbs[1] * dt_size], ncells, MPI_UNSIGNED_CHAR, rank - 1, rank, MPI_COMM_WORLD, &reqs[0]);
        MPI_Irecv(&buf[0], ncells, MPI_UNSIGNED_CHAR, rank - 1, rank - 1, MPI_COMM_WORLD, &reqs[1]);
        
        MPI_Wait(&reqs[0], &stats[0]);
        MPI_Wait(&reqs[1], &stats[1]);
    }
    
    else
    {
        // For others, they exchange both
        MPI_Irecv(&buf[0], ncells, MPI_UNSIGNED_CHAR, rank - 1 , rank - 1, MPI_COMM_WORLD, &reqs[0]);
        MPI_Isend(&buf[border * gbs[0] * gbs[1] * dt_size], ncells, MPI_UNSIGNED_CHAR, rank - 1, rank, MPI_COMM_WORLD, &reqs[1]);
        
        MPI_Isend(&buf[bufsize * dt_size], ncells, MPI_UNSIGNED_CHAR, rank + 1, rank, MPI_COMM_WORLD, &reqs[2]);
        MPI_Irecv(&buf[(bufsize + border * gbs[0] * gbs[1]) * dt_size], ncells, MPI_UNSIGNED_CHAR, rank + 1, rank + 1, MPI_COMM_WORLD, &reqs[3]);
        
        MPI_Waitall(4, reqs, stats);
    }
    
    double communication_end = MPI_Wtime();
    
    
    
 
    
    
    // PERFORM THE ACTUAL BLURING OPERATION
    double compute_start = MPI_Wtime();
    
    // A buffer for the blured image
    unsigned char *bbuffer;
    
    // The number of slices for each process.
    int slices;
    
    bbuffer = malloc(n);
    
    // Take the overflow into consideration
    if(rank == nprocs - 1)
    {
        slices = gbs[2]/nprocs + overflow;
    }
    else
    {
        slices = gbs[2]/nprocs;
    }
    
    
    // Spatial Filtering
    // Skip the border pixels for bluring
    //uint16
    
    /*
    if (dt_size == 2)
    {
        short x = 0;
    
        // Low-passing filtering
        int i = 1, j = 1, k = 1;
        for(k = border; k < slices - border; k++)
        {
            for(i = border; i < gbs[0] - border; i++)
            {
                for(j = border; j < gbs[1] - border; j++)
                {
                    short y = 0;
                    int *ind = generate_filtering_matrix((k + border), i, j);
                    
                    int m = 0;
                    for(m = 0; m < avg_value; m++)
                    {
                        memcpy(&x, &buf[ind[m] * dt_size], sizeof(short));
                        y += x/avg_value;
                    }
                    //bbuffer[(k * gbs[0] * gbs[1] + i * gbs[1] + j) * dt_size] = y;
                    memcpy(&bbuffer[(k * gbs[0] * gbs[1] + i * gbs[1] + j) * dt_size], &y, sizeof(short));
                }
            }
        }
        
        // For the ghost cells from the last rank (rank + 1)
        // The first border slices
        if(rank != 0)
        {
            int k = 0, i = 1, j = 1;
            for(k = 0; k < border; k++)
            {
                for(i = border; i < gbs[0] - border; i++)
                {
                    for(j = border; j < gbs[1] - border; j++)
                    {
                        short y = 0;
                        int *ind = generate_filtering_matrix((k + border), i, j);
                        
                        int m = 0;
                        for(m = 0; m < avg_value; m++)
                        {
                            //bbuffer[i * gbs[1] + j] += buf[ind[m]]/avg_value;
                            memcpy(&x, &buf[ind[m] * dt_size], sizeof(short));
                            y += x/avg_value;
                        }
                        //bbuffer[(k * gbs[0] * gbs[1] + i * gbs[1] + j) * dt_size] = y;
                        memcpy(&bbuffer[(k * gbs[0] * gbs[1] + i * gbs[1] + j) * dt_size], &y, sizeof(short));
                    }
                }
            }
        }
        
        
        
        // For the ghost cells from next rank (rank + 1)
        // The last border slices
        if(rank != nprocs - 1)
        {
            int k = 0, i = 1, j = 1;
            for(k = slices - border; k < slices; k++)
            {
                for(i = border; i < gbs[0] - border; i++)
                {
                    for(j = border; j < gbs[1] - border; j++)
                    {
                        short y = 0;
                        int *ind = generate_filtering_matrix((k + border), i, j);
                        
                        int m = 0;
                        for(m = 0; m < avg_value; m++)
                        {
                            //bbuffer[last_slice + i * gbs[1] + j] += buf[ind[m]]/avg_value;
                            memcpy(&x, &buf[ind[m] * dt_size], sizeof(short));
                            y += x/avg_value;
                        }
                        //bbuffer[(k * gbs[0] * gbs[1] + i * gbs[1] + j) * dt_size] = y;
                        memcpy(&bbuffer[(k * gbs[0] * gbs[1] + i * gbs[1] + j) * dt_size], &y, sizeof(short));
                    }
                }
            }
        }
    }
    
    
    //uint8
    else
    {
     
     */
    
    

    
    /*
    if(dt_size == 2)
    {
        test.s = 0;
    }
    
    printf("The size of bar.value.s is %d \n", sizeof(test.s));
     */
    
    //data1.s = 0;
    
    int i = 1, j = 1, k = 1;
    for(k = border; k < slices - border; k++)
    {
        for(i = border; i < gbs[0] - border; i++)
        {
            for(j = border; j < gbs[1] - border; j++)
            {
                
                unsigned char * ret = generate_filtering_matrix((k + border), i, j, dt_size, buf);
                memcpy(&bbuffer[(k * gbs[0] * gbs[1] + i * gbs[1] + j) * dt_size], ret, dt_size);
                
                 //memcpy(&blur_buf[(k * gbs[0] * gbs[1] + i * gbs[1] + j) * dt_size], &data2.s, dt_size);
                /*
                data2.s = 0;
                int *ind = generate_filtering_matrix((k + border), i, j);
                
                int m = 0;
                for(m = 0; m < avg_value; m++)
                {
                 
                    //memcpy(&bar.value.s, &buf[m * dt_size], dt_size);
                    memcpy(&data1.s, &buf[ind[m] * dt_size], dt_size);
                    data2.s += data1.s/avg_value;
                }
                //bbuffer[(k * gbs[0] * gbs[1] + i * gbs[1] + j * dt_size)] += buf[m]/avg_value;
                memcpy(&bbuffer[(k * gbs[0] * gbs[1] + i * gbs[1] + j) * dt_size], &data2.s, dt_size);
                 */
            }
        }
    }

    // For the ghost cells from last rank (rank -1), which is stored in the first slice.
    // k = 1
    if(rank != 0)
    {
        //int k = 0, i = 1, j = 1;
        for(k = 0; k < border; k++)
        {
            for(i = border; i < gbs[0] - border; i++)
            {
                for(j = border; j < gbs[1] - border; j++)
                {
                    
                    unsigned char * ret = generate_filtering_matrix((k + border), i, j, dt_size, buf);
                    
                    if (ret[0] != 0 )
                        printf ("The value ret[0] is %d of (i %d, j %d, k %d) from rank %d\n", ret[0], i, j, k, rank);
                    
                    memcpy(&bbuffer[(k * gbs[0] * gbs[1] + i * gbs[1] + j) * dt_size], ret, dt_size);
                    /*
                    data2.s = 0;
                    int *ind = generate_filtering_matrix((k + border), i, j);
                    
                    int m = 0;
                    for(m = 0; m < avg_value; m++)
                    {
                         //memcpy(&bar.value.s, &buf[m * dt_size], dt_size);
                        memcpy(&data1.s, &buf[ind[m] * dt_size], dt_size);
                        data2.s += data1.s/avg_value;
                        //bbuffer[(k * gbs[0] * gbs[1] + i * gbs[1] + j * dt_size)] += buf[m]/avg_value;
                    }
                    memcpy(&bbuffer[(k * gbs[0] * gbs[1] + i * gbs[1] + j) * dt_size], &data2.s, dt_size);
                     */
                }
            }
        }
    }



    // For the ghost cells from next rank (rank + 1)
    if(rank != nprocs - 1)
    {
        //int k = 0, i = 1, j = 1;
        for(k = slices - border; k < slices; k++)
        {
            for(i = border; i < gbs[0] - border; i++)
            {
                for(j = border; j < gbs[1] - border; j++)
                {
                    unsigned char * ret = generate_filtering_matrix((k + border), i, j, dt_size, buf);
                    
                    if (ret[0] != 0 )
                        printf ("The value ret[0] is %d of (i %d, j %d, k %d) from rank %d\n", ret[0], i, j, k, rank);
                    
                    memcpy(&bbuffer[(k * gbs[0] * gbs[1] + i * gbs[1] + j) * dt_size], ret, dt_size);
                    /*
                    data2.s = 0;
                    int *ind = generate_filtering_matrix((k + border), i, j);
//
                    int m = 0;
                    for(m = 0; m < avg_value; m++)
                    {
                        memcpy(&data1.s, &buf[ind[m] * dt_size], dt_size);
                        data2.s += data1.s/avg_value;
                         //memcpy(&bar.value.s, &buf[m * dt_size], dt_size);
                         //bbuffer[(k * gbs[0] * gbs[1] + i * gbs[1] + j * dt_size)] += buf[m]/avg_value;
                    }
                    memcpy(&bbuffer[(k * gbs[0] * gbs[1] + i * gbs[1] + j) * dt_size], &data2.s, dt_size);
                     */
                }
            }
        }
    }
       
   // }
    
    double compute_end = MPI_Wtime();
    

    // WRITE THE FILE IN PARALLEL (EXACT OPPOSITE of THE FIRST STEP)
    double io_write_start = MPI_Wtime();
   
    // write_file_name is "bluring_" + read_file_name
    sprintf(write_file_name, "bluring_%s", file_name);
    
    // Open the write file. If this file doesn't exist, then create it
    MPI_File_open(MPI_COMM_WORLD, write_file_name, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh[1]);
    
    // Set the read location for each rank
    MPI_File_seek(fh[1], rank * bufsize * dt_size, MPI_SEEK_SET);
    
    // Write the data in bbuffer into file
    MPI_File_write_all(fh[1], bbuffer, n, MPI_UNSIGNED_CHAR, &status[1]);

    MPI_File_close(&fh[1]);
    
    double io_write_end = MPI_Wtime();
    
    
    free(buf);
    free(bbuffer);
    
    
    // Calculating the total time and other decomposistion time
    double total_time = io_write_end - io_read_start;
    double max_time;
    MPI_Allreduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (max_time == total_time)
    {
        printf("Time take to blur %d x %d x %d image is %f \n Time decomposistion (IO + COMM + COMP + IO): %f + %f + %f + %f \n", gbs[0], gbs[1], gbs[2], max_time, (io_read_end - io_read_start),
               (communication_end - communication_start),
               (compute_end - compute_start),
               (io_write_end - io_write_start));
    }
    

    MPI_Finalize();
    
}



// Parse the arguments
// There are four arguments: -g, -f, -i, -m
static void parse_args(int argc, char **argv)
{
    char flags[] = "g:f:i:m:";
    int one_opt = 0;
    
    while((one_opt = getopt(argc, argv, flags)) != EOF)
    {
        switch (one_opt) {
            case('g'): //global dimention, e.g., 256x256x256
                if((sscanf(optarg, "%dx%dx%d", &gbs[0], &gbs[1], &gbs[2]) == EOF))
                    exit(-1);
                break;
                
            case('f'): // read file name
                if (sprintf(file_name, "%s", optarg) < 0)
                    exit(-1);
                break;
            
            case('i'): // image formate, e.g., uint8
                if((sscanf(optarg, "uint%d", &image_format) == EOF) || image_format > 64)
                    exit(-1);
                break;
                
            case('m'): // the size of matrix, e.g., 3 means 3 * 3 * 3 matrix
                if((sscanf(optarg, "%d", &matrix_size) == EOF || matrix_size % 2 == 0))
                    exit(-1);
                break;
        }
    }
}

