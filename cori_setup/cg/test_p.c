/*
*   Matrix Market I/O example program
*
*   Read a real (non-complex) sparse matrix from a Matrix Market (v. 2.0) file.
*   and copies it to stdout.  This porgram does nothing useful, but
*   illustrates common usage of the Matrix Matrix I/O routines.
*   (See http://math.nist.gov/MatrixMarket for details.)
*
*   Usage:  a.out [filename] > output
*
*
*   NOTES:
*
*   1) Matrix Market files are always 1-based, i.e. the index of the first
*      element of a matrix is (1,1), not (0,0) as in C.  ADJUST THESE
*      OFFSETS ACCORDINGLY offsets accordingly when reading and writing
*      to files.
*
*   2) ANSI C requires one to use the "l" format modifier when reading
*      double precision floating point numbers in scanf() and
*      its variants.  For example, use "%lf", "%lg", or "%le"
*      when reading doubles, otherwise errors will occur.
*/

#include <string.h>
#include <malloc.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"
#include "poski.h"
#include <time.h>
#include "utils.h"
#include "poski_partitioncommon.h"
#include "poski_kernelcommon.h"
#include "poski_matrixcommon.h"
#include "poski_vectorcommon.h"
#include "poski_malloc.h"
#include "poski_print.h"


double dot_prod(double* a,double* b,int n){
   double ret =0.0;
   //#pragma omp for
   for (int i=0;i<n;i++){
      //printf("residue %d %lf %lf %lf: \n",i,ret,a[i],b[i]);
      ret+= a[i]*b[i];
   }
   return ret;
}

// target = target +/- mult * to_process
void vec_add(double* target, double* source ,double* to_process, double mult, int n, int add)
{
   mult = mult*add;
   //#pragma omp for
   for (int i=0;i<n;i++){
      target[i] = source[i] + mult*to_process[i];
   }   
}


int main(int argc, char *argv[])
{
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;
    int i;
    srand((unsigned) 1);
    if (argc < 2)
	{
		fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
		exit(1);
	}
    else
    {
        if ((f = fopen(argv[1], "r")) == NULL)
            exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    int n_threads = atoi(argv[2]);

    //scanf("%d",argv[2]); 
    //printf("%d\n",n_threads); 
    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        exit(1);


    /* reseve memory for matrices */

    int* II = (int *) malloc(nz * sizeof(int));
    int* JJ = (int *) malloc(nz * sizeof(int));
    int* JA= (int *) malloc((M+1) * sizeof(int));
    double* val = (double *) malloc(nz * sizeof(double));


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */
    JA[0]=0;
    int JA_idx=1;
    for (i=0; i<nz; i++)
    {
        fscanf(f, "%d %d %lg\n", &II[i], &JJ[i], &val[i]);
        II[i]--;  /* adjust from 1-based to 0-based */
        JJ[i]--;
        if (II[i]==JJ[i]){
           val[i] = val[i]/2;
           // since the MM format only stores half of the matrix when symmetric, we divide the diagnoal by half, and reulting M satisfies A = M + M.T
        }
        //if (i>0 && J[i]>J[i-1]){
        //   JA_idx+=1;
        //   JA[JA_idx]=i;
        //}
        while (JJ[i]>=JA_idx){
            JA[JA_idx]=i;
            JA_idx++;
        }
    }
    JA[M] = nz;

    if (f !=stdin) fclose(f);

    /************************/
    /* now write out matrix */
    /************************/

    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);
    //for (i=0; i<nz; i++)
    //    fprintf(stdout, "%d %d %20.19g\n", I[i], J[i], val[i]);
    time_t s;
    struct tm* current_time;

    printf("begin pOski!\n");
    s = time(NULL);   
    // to get current time 
    current_time = localtime(&s); 
    // print time in minutes, 
    // hours and seconds 
    printf("start at %02d:%02d:%02d\n", current_time->tm_hour, current_time->tm_min, current_time->tm_sec);

        //STEP 1: /* Initialize pOSKI library */
        poski_Init();
        //STEP 2: /* Initialize Sparse matrix A in CSR format, and dense vectors */
        int nrows=M; int ncols=N; int nnz=nz;
        //int Aptr[4]={0, 1, 3, 5}; int Aind[5]={0, 0, 1, 0, 2}; double Aval[5]={1, -2, 1, 0.5, 1};
        double *x = (double *) malloc(M * sizeof(double));
        double *r = (double *) malloc(M * sizeof(double));
        double *p = (double *) malloc(M * sizeof(double));
        double *b = (double *) malloc(M * sizeof(double));
        double *sol = (double *) malloc(M * sizeof(double));
        double *Ap = (double *) malloc(M * sizeof(double));


        //STEP 3: /* Create a default thread object {with #threads = #available_cores} */
        clock_t before = clock();
        poski_threadarg_t *poski_thread = poski_InitThreads();
        //poski_partitionarg_t *partitionMat =1;
        //poski_partitionarg_t *partitionMat = poski_PartitionMatHints(SemiOneD, n_threads, poski_MatMult,OP_NORMAL);
        //printf("%p",poski_PartitionMatHints(SemiOneD, n_threads, NULL,OP_NORMAL));

        poski_ThreadHints(poski_thread, NULL, POSKI_PTHREAD, n_threads);
        //poski_sparse_matrix_t *Sp_A = poski_LoadMatrix_MM_to_CSR (, poski_index_t MakeUnsymmetric)


        //STEP 4: /* Create a tunable-matrix object by wrapping the partitioned sub-matrices using a thread object and a default partition-matrix object {with #partitions = #threads} */
        poski_mat_t A_tunable =
        poski_CreateMatCSR ( JA, II, val, nrows, ncols, nnz,/* Sparse matrix A in CSR format */
        SHARE_INPUTMAT, /* <matrix copy mode> */
        poski_thread, /* <thread object> */
        NULL,//partitionMat, /* <partition-matrix object> (NULL: default) */
        2, INDEX_ZERO_BASED, MAT_GENERAL);/* specify how to interpret non-zero pattern */



        //poski_ThreadHints(NULL, A_tunable, POSKI_OPENMP, n_threads);

        //STEP 5: /* Create wrappers around the dense vectors with <partition-vector object> (NULL: default) */
        poski_vec_t x_view = poski_CreateVec(x, nrows, STRIDE_UNIT, NULL);
        poski_vec_t r_view = poski_CreateVec(r, nrows, STRIDE_UNIT, NULL);
        poski_vec_t p_view = poski_CreateVec(p, nrows, STRIDE_UNIT, NULL);
        poski_vec_t b_view = poski_CreateVec(b, nrows, STRIDE_UNIT, NULL);
        poski_vec_t sol_view = poski_CreateVec(sol, nrows, STRIDE_UNIT, NULL);
        poski_vec_t Ap_view = poski_CreateVec(Ap, nrows, STRIDE_UNIT, NULL);

        //STEP 6: /* Partition input/output vectors and Perform matrix vector multiply (SpMV), y = α!Ax + β!y */
        //poski_report_MatMultCSR(A_tunable, OP_NORMAL, alpha, x_view, beta, y_view);
        //for (int temp=0;temp<100000;temp++)poski_MatMult(A_tunable, OP_NORMAL, alpha, x_view, beta, y_view);
        

        s = time(NULL);
        // to get current time 
        current_time = localtime(&s);
        // print time in minutes, 
        // hours and seconds 
        printf("finish init at %02d:%02d:%02d\n", current_time->tm_hour, current_time->tm_min, current_time->tm_sec);


        // begin actual algorithm
        //initialization find sol,b
        //#pragma omp for
        for (int temp=0;temp<nrows;temp++) {
            sol[temp]=temp;
            x[temp]=0;r[temp]=0;b[temp]=0.0;
        }
        poski_MatMult(A_tunable, OP_NORMAL, 1.0, sol_view, 1.0, b_view);
        poski_MatMult(A_tunable, OP_TRANS, 1.0, sol_view, 1.0, b_view);
        
        //#pragma omp for
        for (int temp=0;temp<nrows;temp++) {
            r[temp]=b[temp];
            p[temp]=b[temp];
        }
        // r = b - A*x, since x is 0 initially, p=r=b
        
        double residue = dot_prod(r,r,M),prev_residue;
        double alpha,beta;
        // residue = r.T @ r
        residue = dot_prod(r,r,M);printf("current residue is %lf. \n",residue);

        current_time = localtime(&s);
        printf("begin while loop at %02d:%02d:%02d\n", current_time->tm_hour, current_time->tm_min, current_time->tm_sec);
        int it=0;
        while (it<10000){
            it++;
            // Ap = A @ p
            poski_MatMult(A_tunable, OP_NORMAL, 1.0, p_view, 0.0, Ap_view);
            poski_MatMult(A_tunable, OP_TRANS, 1.0, p_view, 1.0, Ap_view);

           /* alpha = residue / dot_prod(p,Ap,M) ;
            //printf(" current alpha is %lf\n",alpha);

            // x = x + alpha*p
            vec_add(x,x,p,alpha,M,1);
*/

            if ((it+1)%20==0){
               //printf("current residue is %lf. \n",residue);
               //r = b-A@x
               for (int temp=0;temp<M;temp++)r[temp]=b[temp];
               poski_MatMult(A_tunable, OP_NORMAL,-1.0, x_view, 1.0, r_view);
               poski_MatMult(A_tunable, OP_TRANS, -1.0, x_view, 1.0, r_view);
            }
            else{
                // r = r - alpha* Ap
                vec_add(r,r,Ap,alpha,M,-1); 
            }
            prev_residue = residue;
            // residue = r.T @ r
            residue = dot_prod(r,r,M);
            
            beta = residue / prev_residue; 
            //printf(" current beta is %lf\n",beta);

            vec_add(p,r,p,beta,M,1);
        }




        clock_t difference = clock() - before;
        s = time(NULL);
        // to get current time 
        current_time = localtime(&s);
        // print time in minutes, 
        // hours and seconds 
        printf("end at %02d:%02d:%02d\n", current_time->tm_hour, current_time->tm_min, current_time->tm_sec);
        printf("it takes %f seconds to finish!\n\n",1.0*difference / CLOCKS_PER_SEC);

        //for (int temp=0;temp<nrows;temp++)printf("%d %f\n",temp,*(y_view->vec+temp));

        //STEP 7: /* Clean-up interface objects and threads, and shut down pOSKI library */
        poski_DestroyMat(A_tunable); 
        //poski_DestroyVec(x_view);
        //poski_DestroyVec(r_view);
        //poski_DestroyVec(p_view);
        //poski_DestroyVec(b_view);
        //poski_DestroyVec(sol_view);
        //poski_DestroyVec(Ap_view);
        //poski_DestroyThreads(poski_thread);
        poski_Close();
	return 0;
}
