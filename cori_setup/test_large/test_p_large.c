#include <stdio.h>
#include "poski.h"
int main(int argc, char **argv)
{
	//STEP 1: /* Initialize pOSKI library */
	poski_Init();
	//STEP 2: /* Initialize Sparse matrix A in CSR format, and dense vectors */
	int nrows=3; int ncols=3; int nnz=5;
	int Aptr[4]={0, 1, 3, 5}; int Aind[5]={0, 0, 1, 0, 2}; double Aval[5]={1, -2, 1, 0.5, 1};
	double x[3]={.25, .45, .65}; double y[3]={1, 1, 1};
	double alpha = -1, beta = 1;

	//STEP 3: /* Create a default thread object {with #threads = #available_cores} */
	poski_threadarg_t *poski_thread = poski_InitThreads();

	//STEP 4: /* Create a tunable-matrix object by wrapping the partitioned sub-matrices using a thread object and a default partition-matrix object {with #partitions = #threads} */
	poski_mat_t A_tunable = 
	poski_CreateMatCSR ( Aptr, Aind, Aval, nrows, ncols, nnz,/* Sparse matrix A in CSR format */
	SHARE_INPUTMAT, /* <matrix copy mode> */
	poski_thread, /* <thread object> */
	NULL, /* <partition-matrix object> (NULL: default) */
	2, INDEX_ZERO_BASED, MAT_GENERAL);/* specify how to interpret non-zero pattern */


	//STEP 5: /* Create wrappers around the dense vectors with <partition-vector object> (NULL: default) */
	poski_vec_t x_view = poski_CreateVec(x, 3, STRIDE_UNIT, NULL);
	poski_vec_t y_view = poski_CreateVec(y, 3, STRIDE_UNIT, NULL);

	//STEP 6: /* Partition input/output vectors and Perform matrix vector multiply (SpMV), y = α!Ax + β!y */
	//poski_report_MatMultCSR(A_tunable, OP_NORMAL, alpha, x_view, beta, y_view);
	printf("%d\n",poski_MatMult(A_tunable, OP_NORMAL, alpha, x_view, beta, y_view));
	printf("%f %f %f\n",*y_view->vec,*(y_view->vec+1),*(y_view->vec+2));

	//STEP 7: /* Clean-up interface objects and threads, and shut down pOSKI library */
	poski_DestroyMat(A_tunable); poski_DestroyVec(x_view); poski_DestroyVec(y_view);
	poski_DestroyThreads(poski_thread);
	poski_Close();

	return 0;
}
