/**
 * C program to solve the two-dimensional Poisson equation on
 * a unit square using one-dimensional eigenvalue decompositions
 * and fast sine transforms.
 *
 * Einar M. Rønquist
 * NTNU, October 2000
 * Revised, October 2001
 * Revised by Eivind Fonn, February 2015
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#define PI 3.14159265358979323846
#define true 1
#define false 0

typedef double real;
typedef int bool;

// Function prototypes
real *mk_1D_array(size_t n, bool zero);
real **mk_2D_array(size_t n1, size_t n2, bool zero);
void transpose(real **bt, real **b, size_t m);
real rhs(real x, real y);
void print_matrix(double **b, int m, int n);
void parallel_transpose(real **bt, real **b, real *send, real *recv, int world_size, int m, int n, int world_rank);
int is_pow2(int n);

// Functions implemented in FORTRAN in fst.f and called from C.
// The trailing underscore comes from a convention for symbol names, called name
// mangling: if can differ with compilers.
void fst_(real *v, int *n, real *w, int *nn);
void fstinv_(real *v, int *n, real *w, int *nn);

int main(int argc, char **argv)
{
    if (argc < 3) {
        printf("Not enough arguments:\n");
        printf("./poisson <problem size> <number of threads per process>\n");
    }

    int t = atoi(argv[2]);
    omp_set_num_threads(t);



    int world_size, world_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    //printf("World size: %i\n", world_size);

    //Check if world_size is a power of 2
    printf("w_size: %i\n", world_size);
    //real l = log(world_size)/log(2);

    /*
     *  The equation is solved on a 2D structured grid and homogeneous Dirichlet
     *  conditions are applied on the boundary:
     *  - the number of grid points in each direction is n+1,
     *  - the number of degrees of freedom in each direction is m = n-1,
     *  - the mesh size is constant h = 1/n.
     */
    int n = atoi(argv[1]);
    int m = n/world_size;
    real h = 1.0 / n;


  	if (!is_pow2(world_size)) {
  		if (world_rank == 0) {
  			printf("Number of processes is not a power of 2. Aborting\n");
  		}
  		MPI_Abort(MPI_COMM_WORLD, MPI_ERR_DIMS);
  		return 0;
  	} else if (!is_pow2(n)) {
      if (world_rank == 0) {
  			printf("Problem size is not a power of 2. Aborting\n");
  		}
  		MPI_Abort(MPI_COMM_WORLD, MPI_ERR_DIMS);
  		return 0;
    }

    real e_time;
    if (world_rank == 0) {
      e_time = MPI_Wtime();
    }



    real **b = mk_2D_array(m, n, false);
    real *send = mk_1D_array((size_t)(m*n), false);
    real *recv = mk_1D_array((size_t)(m*n), false);
    real **bt = mk_2D_array(m, n, false);

    int nn = 4 * n;
    real *z = mk_1D_array(nn, false);

    /*
     * Grid points are generated with constant mesh size on both x- and y-axis.
     */
    real *x_grid = mk_1D_array(m, false);
    real *y_grid = mk_1D_array(n+1, false);

    //#pragma omp parallel for
    for (size_t i = 0; i < m; i++) {
        x_grid[i] = (i + m*world_rank) * h ;
    }

    /*if (world_rank == 3) {
      for (int i=0; i<m;i++) {
        printf("%lf\n", x_grid[i]);
      }
    }*/

    //#pragma omp parallel for
    for (size_t i = 0; i < n+1; i++) {
        y_grid[i] = i * h ;
    }

    /*if (world_rank == 0) {
      for (int i=0; i<n+1;i++) {
        printf("%lf\n", y_grid[i]);
      }
    }*/

    real *diag = mk_1D_array(n, false);
    //#pragma omp parallel for
    for (size_t i = 0; i < n; i++) {
        diag[i] = 2.0 * (1.0 - cos((i+1) * PI / n));
    }
    int debug_rank = 0;

    //if (world_rank == debug_rank) printf("diag: %lf %lf %lf %lf \n", diag[0], diag[1], diag[2], diag[3]);


    //#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            b[i][j] = h * h * rhs(x_grid[i+1], y_grid[j+1]);
            //b[i][j] = rhs(x_grid[i+1], y_grid[j+1]); for easier testing
        }
    }

    /*if (world_rank == 0) {
      //printf("%i\n", world_rank);
      print_matrix(b, m, n);
    }*/


    /*
     * Compute \tilde G^T = S^-1 * (S * G)^T (Chapter 9. page 101 step 1)
     * Instead of using two matrix-matrix products the Discrete Sine Transform
     * (DST) is used.
     * The DST code is implemented in FORTRAN in fsf.f and can be called from C.
     * The array zz is used as storage for DST coefficients and internally for
     * FFT coefficients in fst_ and fstinv_.
     * In functions fst_ and fst_inv_ coefficients are written back to the input
     * array (first argument) so that the initial values are overwritten.
     */



    //#pragma omp parallel for
    for (size_t i = 0; i < m; i++) {
        fst_(b[i], &n, z, &nn);
    }

    /*if (world_rank == debug_rank) {
      //printf("%i\n", world_rank);
      print_matrix(b, m, n-1);
    }*/

    parallel_transpose(bt, b, send, recv, world_size, m, n, world_rank);

    /*if (world_rank == debug_rank) {
      printf("\n");
      print_matrix(bt, m, n-1);
    }*/


    //#pragma omp parallel for
    for (size_t i = 0; i < m; i++) {
        fstinv_(bt[i], &n, z, &nn);
    }

    /*if (world_rank == debug_rank) {
      printf("\n");
      print_matrix(bt, m, n-1);
    }*/


    /*
     * Solve Lambda * \tilde U = \tilde G (Chapter 9. page 101 step 2)
     */

    //#pragma omp parallel for collapse(2)
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            bt[i][j] = bt[i][j] / (diag[i+world_rank*m] + diag[j]);
        }
    }

    /*if (world_rank == debug_rank) {
      printf("\n");
      print_matrix(bt, m, n-1);
    }*/

    /*
     * Compute U = S^-1 * (S * Utilde^T) (Chapter 9. page 101 step 3)
     */

    //#pragma omp parallel for
    for (size_t i = 0; i < m; i++) {
        fst_(bt[i], &n, z, &nn);
    }
    parallel_transpose(b, bt, send, recv, world_size, m, n, world_rank);
    //transpose(b, bt, m);


    //#pragma omp parallel for
    for (size_t i = 0; i < m; i++) {
        fstinv_(b[i], &n, z, &nn);
    }

    /*
     * Compute maximal value of solution for convergence analysis in L_\infty
     * norm.
     */

    //print_matrix(b, m, m);




      //if(world_rank != world_size-1) {
        double u_max = 0.0;
        //size_t o = world_rank == 0 ? 1 : 0;
        //#pragma omp parallel for collapse(2)
        //for (size_t i = world_rank == 0 ? 1 : 0; i < m; i++) {
        for (size_t i = 0; i < m; i++) {
            for (size_t j = 0; j < n-1; j++) {
                u_max = u_max > b[i][j] ? u_max : b[i][j];
            }
        }

        //printf("u_max = %e\n", u_max);

        real u_max_glob;
        MPI_Reduce(&u_max, &u_max_glob, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        if (world_rank == 0) printf("u_max = %e\n", u_max_glob);
    //}

    MPI_Finalize();
    return 0;
}

/*
 * This function is used for initializing the right-hand side of the equation.
 * Other functions can be defined to swtich between problem definitions.
 */

real rhs(real x, real y) {
    return 2 * (y - y*y + x - x*x);
    //return x;
    //return 1;
}


void transpose(real **bt, real **b, size_t m)
{
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            bt[i][j] = b[j][i];
        }
    }
}

/*
 * The allocation of a vectore of size n is done with just allocating an array.
 * The only thing to notice here is the use of calloc to zero the array.
 */

real *mk_1D_array(size_t n, bool zero)
{
    if (zero) {
        return (real *)calloc(n, sizeof(real));
    }
    return (real *)malloc(n * sizeof(real));
}

/*
 * The allocation of the two-dimensional array used for storing matrices is done
 * in the following way for a matrix in R^(n1*n2):
 * 1. an array of pointers is allocated, one pointer for each row,
 * 2. a 'flat' array of size n1*n2 is allocated to ensure that the memory space
 *   is contigusous,
 * 3. pointers are set for each row to the address of first element.
 */

real **mk_2D_array(size_t n1, size_t n2, bool zero)
{
    // 1
    real **ret = (real **)malloc(n1 * sizeof(real *));

    // 2
    if (zero) {
        ret[0] = (real *)calloc(n1 * n2, sizeof(real));
    }
    else {
        ret[0] = (real *)malloc(n1 * n2 * sizeof(real));
    }

    // 3
    for (size_t i = 1; i < n1; i++) {
        ret[i] = ret[i-1] + n2;
    }
    return ret;
}

void print_matrix(real **b, int m, int n) {
   for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
         printf("%f ", b[i][j]);
      }
      printf("\n");
   }
}

void parallel_transpose(real **bt, real **b, real *send, real *recv, int world_size, int m, int n, int world_rank) {
  size_t i, j;
  #pragma omp parallel for collapse(2)
  for (i=0; i<(size_t)m; i++) {
    for (j=0; j<(size_t)n; j++) {
      send[m*i + (j/m)*(m*m) + j%m] = b[i][j];
      /*if (world_rank == 1) {
        printf("%i\n", (int)(world_size*i + j/world_size*(world_size*world_size) + j%world_size));
      }*/
    }
  }

  //printf("Send, %i: %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", world_rank, send[0], send[1], send[2], send[3], send[4], send[5], send[6], send[7], send[8], send[9], send[10], send[11], send[12], send[13], send[14], send[15]);

  //MPI_Alltoall(&send, m*m, MPI_DOUBLE, &recv, m*m, MPI_DOUBLE, MPI_COMM_WORLD);
  //MPI_Barrier(MPI_COMM_WORLD);
  MPI_Alltoall(&send[0], m*m, MPI_DOUBLE, &recv[0], m*m, MPI_DOUBLE, MPI_COMM_WORLD);

  //printf("Recv, %i: %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", world_rank, recv[0], recv[1], recv[2], recv[3], recv[4], recv[5], recv[6], recv[7], recv[8], recv[9], recv[10], recv[11], recv[12], recv[13], recv[14], recv[15]);

  int cnt = 0;
  //#pragma omp parallel for collapse(2)
  for (j=0; j<(size_t)n; j++) {
    for (i=0; i<(size_t)m; i++) {
      bt[i][j] = recv[cnt];
      cnt++;
    }
  }

  /*if (world_rank == 0) {
    print_matrix(bt, m, n);
  }*/


}

int is_pow2(int n) {
  real l = log(n)/log(2);
  return floor(l) == l;
}
