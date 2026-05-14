#include "pch.h"
#include "framework.h"
#include "matrix.h"
#include <stdlib.h>


//const double RADIX = 2;
//const double LUTINY = 1E-20;

static double* cyno_alloc_double_array(long size) {
	if (size <= 0) {
		return NULL;
	}
	return (double*)malloc((size_t)size * sizeof(double));
}

static long* cyno_alloc_long_array(long size) {
	if (size <= 0) {
		return NULL;
	}
	return (long*)malloc((size_t)size * sizeof(long));
}

static double** cyno_alloc_matrix(long rows, long cols) {
	double** matrix;
	double* data;
	long i;

	if (rows <= 0 || cols <= 0) {
		return NULL;
	}

	matrix = (double**)malloc((size_t)rows * sizeof(double*));
	if (matrix == NULL) {
		return NULL;
	}

	data = (double*)malloc((size_t)(rows * cols) * sizeof(double));
	if (data == NULL) {
		free(matrix);
		return NULL;
	}

	for (i = 0; i < rows; i++) {
		matrix[i] = data + (i * cols);
	}

	return matrix;
}

static void cyno_free_matrix(double** matrix) {
	if (matrix != NULL) {
		free(matrix[0]);
		free(matrix);
	}
}

static void cyno_copy_flat_to_matrix(const double* input, long rows, long cols, double** matrix) {
	long i, j;

	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			matrix[i][j] = input[i * cols + j];
		}
	}
}

static void cyno_copy_matrix_to_flat(double** matrix, long rows, long cols, double* output) {
	long i, j;

	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			output[i * cols + j] = matrix[i][j];
		}
	}
}

bool _stdcall cyno_DiagonalMatrixMultiply(long matrixOneDim, double* matrixOne, long matrixTwoDim1, long matrixTwoDim2, double** matrixTwo);
bool _stdcall cyno_DiagonalMatrixMultiplyMatrixVector(long matrixOneDim1, long matrixOneDim2, double** matrixOne, long matrixTwoDim, double* matrixTwo);
bool UnbalanceVectors(long vMatrixDim1, long vMatrixDim2, double** vMatrix, long transDim, double* trans, short bL, short bH, bool isLeft);
bool BalanceMatrix(long theMatrixDim1, long theMatrixDim2, double** theMatrix, long transDim, double* trans, short* bL, short* bH);
bool Hessenberg(long aMatrixDim1, long aMatrixDim2, double** aMatrix, long transDim, double* trans, long zMatrixDim1, long zMatrixDim2, double** zMatrix, short bL, short bH, bool eigvects);
bool HessEigenSystem(long hMatrixDim1, long hMatrixDim2, double** hMatrix, long zMatrixDim1, long zMatrixDim2, double** zMatrix, long wrDim, double* wr, long wiDim, double* wi, short bL, short bH, bool eigvects);

// It's caller's responsibility to resize output matrix properly save processing time
// and avoid potential memory leaking
// Right most dimension in SAFEARRAY is stored first in array element area.
// For example, in VBA
// Dim A(1 to 2, 1 to 3) as double
// The order of elements in SAFEARRAY element area is
// (1,1), (2,1), (1,2), (2, 2), (1,3), (2,3)

//Purpose
//		Multiply a scalar value to a one dimension or two dimension matrix
//Parameters
//		a			input	scalar value
//		b			input	1D/2D matrix
//		c			output	1D/2D result matrix
//Return Values
//		bool
//Remarks
//		Arrays b and c must be sized properly and in same size before call this function.

CYNOMATHUTILITY_API bool _stdcall cyno_MatrixMultiply(double a, SAFEARRAY** b, SAFEARRAY** c) {
	HRESULT lResult; // return code for OLE functions

	double* pArrayElements_b, * pArrayElements_c; // pointer to the elements of the array
	long i, j;

	bool ret;

	long D1Elements, D2Elements;


	pArrayElements_b = (double*)(*b)->pvData;
	pArrayElements_c = (double*)(*c)->pvData;

	D1Elements = (*b)->rgsabound[0].cElements;
	D2Elements = (*b)->rgsabound[1].cElements;


	lResult = SafeArrayLock(*b);
	if (lResult) return false;
	lResult = SafeArrayLock(*c);
	if (lResult) return false;

	if ((*b)->cDims == 1) {

		for (i = 0; i < D1Elements; i++)
			pArrayElements_c[i] = a * pArrayElements_b[i];
		ret = true;
	}
	else if ((*b)->cDims == 2) {
		for (i = 0; i < D1Elements; i++)
			for (j = 0; j < D2Elements; j++) {
				pArrayElements_c[i * D1Elements + j] = a * pArrayElements_b[i * D1Elements + j];
			}
		ret = true;
	}
	else
		ret = false;
	lResult = SafeArrayUnlock(*b);
	if (lResult) return false;
	lResult = SafeArrayUnlock(*c);
	if (lResult) return false;
	return ret;
}

//Purpose
//		Multiply one or two dimension matrix to another one dimension or two dimension matrix
//Parameters
//		a			input	1D/2D matrix
//		b			input	1D/2D matrix
//		c			output	1D/2D result matrix
//Return Values
//		bool
//Remarks
//		Arrays a, b, c must be sized properly. Array a and b must multiplyable.

CYNOMATHUTILITY_API bool _stdcall cyno_MatrixMultiply_2(SAFEARRAY** a, SAFEARRAY** b, SAFEARRAY** c) {
	HRESULT lResult; 
	double* pArrayElements_a, * pArrayElements_b, * pArrayElements_c; // pointer to the elements of the array

	bool ret = false;

	int aD1Elements, aD2Elements;
	int bD1Elements, bD2Elements;

	int D1, D2, i, j, k;
	double tmp;



	D1 = (*a)->cDims;
	D2 = (*b)->cDims;

	pArrayElements_a = (double*)(*a)->pvData;
	pArrayElements_b = (double*)(*b)->pvData;
	pArrayElements_c = (double*)(*c)->pvData;

	aD1Elements = (*a)->rgsabound[0].cElements;
	bD1Elements = (*b)->rgsabound[0].cElements;

	lResult = SafeArrayLock(*a);
	if (lResult) return false;
	lResult = SafeArrayLock(*b);
	if (lResult) return false;
	lResult = SafeArrayLock(*c);
	if (lResult) return false;



	if (D1 == 2 && D2 == 1) {
		aD2Elements = (*a)->rgsabound[1].cElements;

		if (aD1Elements == bD1Elements) {
			for (i = 0; i < aD2Elements; i++) { //row
				tmp = 0;
				for (k = 0; k < aD1Elements; k++) //column
					tmp = tmp + pArrayElements_a[i + k * aD2Elements] * pArrayElements_b[k];
				pArrayElements_c[i] = tmp;
			}
			ret = true;
		}
	}
	else if (D1 == 1 && D2 == 2) {
		bD2Elements = (*b)->rgsabound[1].cElements;
		if (aD1Elements == bD2Elements) {
			for (i = 0; i < bD1Elements; i++) {
				tmp = 0;
				for (k = 0; k < aD1Elements; k++)
					//tmp = tmp + A(k) * b(k, i)
					tmp = tmp + pArrayElements_a[k] * pArrayElements_b[i * bD2Elements + k];
				pArrayElements_c[i] = tmp;
			}
			ret = true;
		}
	}
	else if (D1 == 2 && D2 == 2) {
		aD2Elements = (*a)->rgsabound[1].cElements;
		bD2Elements = (*b)->rgsabound[1].cElements;
		if (aD1Elements == bD2Elements) {
			for (i = 0; i < aD2Elements; i++)
				for (j = 0; j < bD1Elements; j++) {
					tmp = 0;
					for (k = 0; k < aD1Elements; k++) {
						//tmp = tmp + A(i,k) * b(k,j)
						tmp = tmp + pArrayElements_a[k * aD2Elements + i] * pArrayElements_b[j * bD2Elements + k];
						
					}
					pArrayElements_c[j * aD2Elements + i] = tmp;
				}
			ret = true;
		}
	}

	lResult = SafeArrayUnlock(*a);
	if (lResult) return false;
	lResult = SafeArrayUnlock(*b);
	if (lResult) return false;
	lResult = SafeArrayUnlock(*c);
	if (lResult) return false;

	return ret;
}

//Purpose
//		Multiply one dimension matrix to another one dimension matrix
//Parameters
//		a			input	1D matrix
//		b			input	1D matrix
//		c			output	scalar result
//Return Values
//		bool
//Remarks
//		Arrays a, b must be sized properly and in same size.

CYNOMATHUTILITY_API bool _stdcall cyno_MatrixMultiply_3(SAFEARRAY** a, SAFEARRAY** b, double* c) {
	HRESULT lResult; 
	double* pArrayElements_a, * pArrayElements_b; // pointer to the elements of the array

	bool ret = false;

	long aD1Elements;
	long bD1Elements;

	long D1, D2, i;
	double tmp;


	D1 = (*a)->cDims;
	D2 = (*b)->cDims;

	pArrayElements_a = (double*)(*a)->pvData;
	pArrayElements_b = (double*)(*b)->pvData;

	aD1Elements = (*a)->rgsabound[0].cElements;
	bD1Elements = (*b)->rgsabound[0].cElements;

	lResult = SafeArrayLock(*a);
	if (lResult) return false;
	lResult = SafeArrayLock(*b);
	if (lResult) return false;

	if (D1 == 1 && D2 == 1) {
		if (aD1Elements == bD1Elements) {
			tmp = 0;
			for (i = 0; i < aD1Elements; i++)
				tmp = tmp + pArrayElements_a[i] * pArrayElements_b[i];
			*c = tmp;
			ret = true;
		}
	}

	lResult = SafeArrayUnlock(*a);
	if (lResult) return false;
	lResult = SafeArrayUnlock(*b);
	if (lResult) return false;

	return ret;
}

// ==========================
// Internal matrix functions
// ==========================
bool cyno_ConvertFromSafeArray1D(SAFEARRAY** a, long dim, double* pa) {
	HRESULT lResult; // return code for OLE functions
	double* pArrayElements_a; // pointer to the elements of the array

	long i;

	pArrayElements_a = (double*)(*a)->pvData;

	lResult = SafeArrayLock(*a);
	if (lResult) return false;

	for (i = 0; i < dim; i++)
		pa[i] = pArrayElements_a[i];

	lResult = SafeArrayUnlock(*a);
	if (lResult) return false;

	return true;
}

bool cyno_ConvertFromSafeArray2D(SAFEARRAY** a, long dim1, long dim2, double** pa) {
	HRESULT lResult; // return code for OLE functions
	double* pArrayElements_a; // pointer to the elements of the array

	long i, j;

	pArrayElements_a = (double*)(*a)->pvData;

	lResult = SafeArrayLock(*a);
	if (lResult) return false;

	for (i = 0; i < dim1; i++)
	{
		for (j = 0; j < dim2; j++)
			pa[i][j] = pArrayElements_a[i + j * dim1];
	}

	lResult = SafeArrayUnlock(*a);
	if (lResult) return false;
	return true;
}

bool cyno_ConvertToSafeArray1D(SAFEARRAY** a, long dim, double* pa) {
	HRESULT lResult; // return code for OLE functions
	double* pArrayElements_a; // pointer to the elements of the array

	long i;

	pArrayElements_a = (double*)(*a)->pvData;

	lResult = SafeArrayLock(*a);
	if (lResult) return false;

	for (i = 0; i < dim; i++)
		pArrayElements_a[i] = pa[i];

	lResult = SafeArrayUnlock(*a);
	if (lResult) return false;

	return true;
}

bool cyno_ConvertToSafeArray2D(SAFEARRAY** a, long dim1, long dim2, double** pa) {
	HRESULT lResult; // return code for OLE functions
	double* pArrayElements_a; // pointer to the elements of the array

	long i, j;

	pArrayElements_a = (double*)(*a)->pvData;

	lResult = SafeArrayLock(*a);
	if (lResult) return false;

	for (i = 0; i < dim1; i++)
	{
		for (j = 0; j < dim2; j++)
			pArrayElements_a[i + j * dim1] = pa[i][j];
	}

	lResult = SafeArrayUnlock(*a);
	if (lResult) return false;
	return true;
}


void _stdcall cyno_MatrixMultiScalar1D(double a, long dim, double* b, double* c) {

	long i;

	for (i = 0; i < dim; i++)
		c[i] = a * b[i];
}


void _stdcall cyno_MatrixMultiScalar2D(double a, long dim1, long dim2, double** b, double** c) {
	long i, j;

	for (i = 0; i < dim1; i++)
		for (j = 0; j < dim2; j++)
			c[i][j] = a * b[i][j];
}


bool _stdcall cyno_MatrixMulti2D1D(long adim1, long adim2, double** a, long bdim, double* b, double* c) {
	bool ret = false;

	long i, k;
	double tmp;

	if (adim2 == bdim) {
		for (i = 0; i < adim1; i++) { //row
			tmp = 0;
			for (k = 0; k < adim2; k++) //column
				tmp = tmp + a[i][k] * b[k];
			c[i] = tmp;
		}
		ret = true;
	}
	return ret;
}

bool _stdcall cyno_MatrixMulti2D2D(long adim1, long adim2, double** a, long bdim1, long bdim2, double** b, double** c) {
	bool ret = false;

	long i, j, k;
	double tmp;

	if (adim2 == bdim1) {
		for (i = 0; i < adim1; i++)
			for (j = 0; j < bdim2; j++) {
				tmp = 0;
				for (k = 0; k < adim2; k++) {
					tmp = tmp + a[i][k] * b[k][j];
				}
				c[i][j] = tmp;
			}
		ret = true;
	}
	return ret;
}

bool _stdcall cyno_MatrixMulti1D1D(long adim, double* a, long bdim, double* b, double* c) {
	bool ret = false;

	long i;
	double tmp;

	if (adim == bdim) {
		tmp = 0;
		for (i = 0; i < adim; i++)
			tmp = tmp + a[i] * b[i];
		*c = tmp;
		ret = true;
	}

	return ret;
}




void _stdcall cyno_MatrixTrans(long dim1, long dim2, double** a, double** b) {

	long i, j;

	for (i = 0; i < dim1; i++)
		for (j = 0; j < dim2; j++)
			b[j][i] = a[i][j];
}

void _stdcall cyno_MatrixCopy(long dim1, long dim2, double** a, double** b) {

	long i, j;

	for (i = 0; i < dim1; i++)
		for (j = 0; j < dim2; j++)
			b[i][j] = a[i][j];
}

bool _stdcall cyno_DiagonalMatrixMultiply(long matrixOneDim, double* matrixOne, long matrixTwoDim1, long matrixTwoDim2, double** matrixTwo) {
	long i, j;

	if (matrixOneDim == matrixTwoDim1 && matrixTwoDim1 == matrixTwoDim2) {
		for (i = 0; i < matrixOneDim; i++) {
			for (j = 0; j < matrixOneDim; j++)
				matrixTwo[i][j] = matrixOne[i] * matrixTwo[i][j];
		}
		return true;
	}
	else
		return false;
}

bool _stdcall cyno_DiagonalMatrixMultiplyMatrixVector(long matrixOneDim1, long matrixOneDim2, double** matrixOne, long matrixTwoDim, double* matrixTwo) {
	long i, j;

	if (matrixOneDim2 == matrixTwoDim && matrixOneDim1 == matrixOneDim2) {
		for (i = 0; i < matrixOneDim1; i++) {
			for (j = 0; j < matrixOneDim1; j++)
				matrixOne[i][j] = matrixOne[i][j] * matrixTwo[j];
		}
		return true;
	}
	else
		return false;
}

bool _stdcall cyno_LUDecomp(long aMatrixDim1, long aMatrixDim2, double** aMatrix, long* pVec, short* d) {
	short i, IMax, j, k;
	double BIG, dum, sum, temp;
	double* vv;

	if (aMatrixDim1 == aMatrixDim2) {
		vv = cyno_alloc_double_array(aMatrixDim1);
		if (vv == NULL)
			return false;
		*d = 1;
		for (i = 0; i < aMatrixDim1; i++) {
			BIG = 0;
			for (j = 0; j < aMatrixDim1; j++) {
				temp = fabs(aMatrix[i][j]);
				if (temp > BIG)
					BIG = temp;
			}
			if (BIG == 0) {
				free(vv);
				return false;
			}
			vv[i] = 1 / BIG;
		}
		for (j = 0; j < aMatrixDim1; j++) {
			for (i = 0; i <= j - 1; i++) {
				sum = aMatrix[i][j];
				//for (k = 1; k <= i - 1; k++) - changed from 1 to 0 intentionally to avoid a bug.
				//								This is because C++ starts array subscript from 0 while
				//								the original VBA code starts from 1.
				for (k = 0; k <= i - 1; k++)
					sum = sum - aMatrix[i][k] * aMatrix[k][j];
				aMatrix[i][j] = sum;
			}
			BIG = 0;
			for (i = j; i < aMatrixDim1; i++) {
				sum = aMatrix[i][j];
				//for (k = 1; k <= j - 1; k++)  - changed from 1 to 0 intentionally to avoid a bug.
				//								This is because C++ starts array subscript from 0 while
				//								the original VBA code starts from 1.
				for (k = 0; k <= j - 1; k++)
					sum = sum - aMatrix[i][k] * aMatrix[k][j];
				aMatrix[i][j] = sum;
				dum = vv[i] * fabs(sum);
				if (dum >= BIG) {
					BIG = dum;
					IMax = i;
				}
			}
			if (j != IMax) {
				for (k = 0; k < aMatrixDim1; k++) {
					dum = aMatrix[IMax][k];
					aMatrix[IMax][k] = aMatrix[j][k];
					aMatrix[j][k] = dum;
				}
				*d = -1 * (*d);
				vv[IMax] = vv[j];
			}
			pVec[j] = IMax;
			if (aMatrix[j][j] == 0)
				aMatrix[j][j] = LUTINY;
			if (j != aMatrixDim1 - 1) {
				dum = 1 / aMatrix[j][j];
				for (i = j + 1; i < aMatrixDim1; i++)
					aMatrix[i][j] = aMatrix[i][j] * dum;
			}
		}
		free(vv);
		return true;
	}
	else
		return false;
}

bool _stdcall cyno_LUBackSub(long luMatrixDim1, long luMatrixDim2, double** luMatrix, long* pVec, double* bVec) {
	long i, ii, ip, j;
	double sum;


	if (luMatrixDim1 == luMatrixDim2) {
		//		ii = 0;  - changed from 0 to -1 intentionally to avoid a bug.
						//								This is because C++ starts array subscript from 0 while
						//								the original VBA code starts from 1.
		ii = -1;
		for (i = 0; i < luMatrixDim1; i++) {
			ip = pVec[i];
			sum = bVec[ip];
			bVec[ip] = bVec[i];



			//			if (ii != 0)  - changed from 0 to -1 intentionally to avoid a bug.
							//								This is because C++ starts array subscript from 0 while
							//								the original VBA code starts from 1.
			if (ii != -1)
				for (j = ii; j <= i - 1; j++)
					sum = sum - luMatrix[i][j] * bVec[j];
			else if (sum != 0)
				ii = i;
			bVec[i] = sum;



		}
		for (i = luMatrixDim1 - 1; i >= 0; i--) {
			sum = bVec[i];
			for (j = i + 1; j < luMatrixDim1; j++)
				sum = sum - luMatrix[i][j] * bVec[j];
			bVec[i] = sum / luMatrix[i][i];

		}

		return true;
	}
	else
		return false;
}

CYNOMATHUTILITY_API bool _stdcall cyno_InverseMatrix(long dim1, long dim2, double** sqMatrix, double** invMatrix) 
{
	double* col;
	double** tmpMatrix;
	long* pVec;
	long i, j;
	bool chk;
	short d;


	if (dim1 == dim2) {
		pVec = cyno_alloc_long_array(dim1);
		col = cyno_alloc_double_array(dim1);
		tmpMatrix = cyno_alloc_matrix(dim1, dim1);
		if (pVec == NULL || col == NULL || tmpMatrix == NULL) {
			free(pVec);
			free(col);
			cyno_free_matrix(tmpMatrix);
			return false;
		}
		cyno_MatrixCopy(dim1, dim2, sqMatrix, tmpMatrix);
		chk = cyno_LUDecomp(dim1, dim2, tmpMatrix, pVec, &d);


		if (chk) {
			for (j = 0; j < dim1; j++) {
				for (i = 0; i < dim1; i++)
					col[i] = 0;
				col[j] = 1;
				chk = cyno_LUBackSub(dim1, dim2, tmpMatrix, pVec, col);
				for (i = 0; i < dim1; i++) {
					invMatrix[i][j] = col[i];

				}
			}
			chk = true;
		}
		else
			chk = false;
		cyno_free_matrix(tmpMatrix);
		free(pVec);
		free(col);

		return chk;
	}
	else
		return false;
}

//Purpose
//		Transpose a matrix
//Parameters
//		psaInMatrix	input	2D matrix
//		psaOutMatrix	output	transposed 2D matrix
//Return Values
//		bool
//Remarks
//		Arrays psaInMatrix and psaOutMatrix must be 2 dimension matrix
//		and must be sized properly before call this function.
//
CYNOMATHUTILITY_API bool _stdcall cyno_MatrixTranspose(SAFEARRAY** psaInMatrix, SAFEARRAY** psaOutMatrix) {
	short i;
	long indim1, indim2, outdim1, outdim2;
	double** pInMatrix, ** pOutMatrix;
	bool ret;

	if ((*psaInMatrix)->cDims == 2) {

		indim1 = (*psaInMatrix)->rgsabound[1].cElements;
		indim2 = (*psaInMatrix)->rgsabound[0].cElements;

		outdim1 = (*psaOutMatrix)->rgsabound[1].cElements;
		outdim2 = (*psaOutMatrix)->rgsabound[0].cElements;

		if (indim1 != outdim2 || indim2 != outdim1)
			return false;

		// Allocate memory
		pInMatrix = cyno_alloc_matrix(indim1, indim2);
		pOutMatrix = cyno_alloc_matrix(outdim1, outdim2);
		if (pInMatrix == NULL || pOutMatrix == NULL) {
			cyno_free_matrix(pInMatrix);
			cyno_free_matrix(pOutMatrix);
			return false;
		}

		ret = cyno_ConvertFromSafeArray2D(psaInMatrix, indim1, indim2, pInMatrix);

		cyno_MatrixTrans(indim1, indim2, pInMatrix, pOutMatrix);

		ret = cyno_ConvertToSafeArray2D(psaOutMatrix, outdim1, outdim2, pOutMatrix);

		// De-Allocate memory to prevent memory leak
		cyno_free_matrix(pInMatrix);
		cyno_free_matrix(pOutMatrix);

		return ret;
	}
	return false;
}

//Purpose
//		Matrix Power - matrix^{n}
//Parameters
//		psaInMatrix	input/output	2D square matrix
//		thePower	input	power
//Return Values
//		bool
//Remarks
//		Array psaInMatrix must be 2 dimension square matrix.
//		The result will be returned through array psaInMatrix.

CYNOMATHUTILITY_API bool _stdcall cyno_MatrixPower(SAFEARRAY** psaInMatrix, short thePower) {
	short i, j, rPow;
	long indim1, indim2, outdim1, outdim2;
	double** pInMatrix, ** pTmpMatrix, ** pOutMatrix;
	bool ret;


	ret = false;

	if ((*psaInMatrix)->cDims == 2) {
		indim1 = (*psaInMatrix)->rgsabound[1].cElements;
		indim2 = (*psaInMatrix)->rgsabound[0].cElements;

		outdim1 = indim1;
		outdim2 = indim2;

		if (indim1 != indim2)
			return false;

		// Allocate memory
		pInMatrix = cyno_alloc_matrix(indim1, indim2);
		pTmpMatrix = cyno_alloc_matrix(indim1, indim2);
		pOutMatrix = cyno_alloc_matrix(outdim1, outdim2);
		if (pInMatrix == NULL || pTmpMatrix == NULL || pOutMatrix == NULL) {
			cyno_free_matrix(pInMatrix);
			cyno_free_matrix(pTmpMatrix);
			cyno_free_matrix(pOutMatrix);
			return false;
		}


		ret = cyno_ConvertFromSafeArray2D(psaInMatrix, indim1, indim2, pInMatrix);


		if (ret) {
			// copy InMatrix to TmpMatrix
			for (i = 0; i < indim1; i++)
				for (j = 0; j < indim2; j++)
					pTmpMatrix[i][j] = pInMatrix[i][j];

			rPow = 1;

			ret = true;

			while (ret && rPow < thePower) {
				ret = cyno_MatrixMulti2D2D(indim1, indim2, pTmpMatrix, indim1, indim2, pInMatrix, pOutMatrix);
				rPow++;
				// copy OutMatrix to TmpMatrix
				for (i = 0; i < indim1; i++)
					for (j = 0; j < indim2; j++)
						pTmpMatrix[i][j] = pOutMatrix[i][j];
			}

			if (ret)
				ret = cyno_ConvertToSafeArray2D(psaInMatrix, outdim1, outdim2, pOutMatrix);
		}

		cyno_free_matrix(pInMatrix);
		cyno_free_matrix(pTmpMatrix);
		cyno_free_matrix(pOutMatrix);
	}
	return ret;
}

//Purpose
//		Cholesky Decompostion
//Parameters
//		psaInMatrix	input/output	2D square matrix
//Return Values
//		bool
//Remarks
//		Array psaInMatrix must be 2 dimension square matrix.
//		The result will be returned through array psaInMatrix.

CYNOMATHUTILITY_API bool _stdcall cyno_CholeskyDecomp(SAFEARRAY** psaSquareMatrix) {
	short i, j, k;
	bool chk;
	long indim1, indim2;
	double** pTmp, ** pSquareMatrix, sum;

	chk = false;

	if ((*psaSquareMatrix)->cDims == 2) {
		indim1 = (*psaSquareMatrix)->rgsabound[1].cElements;
		indim2 = (*psaSquareMatrix)->rgsabound[0].cElements;

		if (indim1 == indim2) {
			pTmp = cyno_alloc_matrix(indim1, indim1);
			pSquareMatrix = cyno_alloc_matrix(indim1, indim1);
			if (pTmp == NULL || pSquareMatrix == NULL) {
				cyno_free_matrix(pTmp);
				cyno_free_matrix(pSquareMatrix);
				return false;
			}

			chk = cyno_ConvertFromSafeArray2D(psaSquareMatrix, indim1, indim1, pSquareMatrix);
			for (i = 0; i < indim1; i++)
				for (j = i; j < indim1; j++) {
					sum = pSquareMatrix[i][j];
					for (k = i - 1; k >= 0; k--)
						sum = sum - pSquareMatrix[i][k] * pSquareMatrix[j][k];
					if (i == j) {
						if (sum < 0) {
							//deallocate memory
							cyno_free_matrix(pSquareMatrix);
							cyno_free_matrix(pTmp);
							return false;
						}
						else
							pTmp[i][i] = sqrt(sum);
					}
					else {
						pTmp[j][i] = sum / pTmp[i][i];
						pSquareMatrix[j][i] = pTmp[j][i];
					}
				}
			for (i = 0; i < indim1; i++)
				for (j = i + 1; j < indim1; j++)
					pTmp[i][j] = 0;

			//copy matrix from pTmp to pSquareMatrix
			for (i = 0; i < indim1; i++)
				for (j = 0; j < indim1; j++)
					pSquareMatrix[i][j] = pTmp[i][j];
			chk = cyno_ConvertToSafeArray2D(psaSquareMatrix, indim1, indim1, pSquareMatrix);

			//deallocate memory
			cyno_free_matrix(pSquareMatrix);
			cyno_free_matrix(pTmp);
			chk = true;
		}
	}

	return chk;
}

double Pythag(double xArg, double yArg)
{
	double absX, absY;
	absX = fabs(xArg);
	absY = fabs(yArg);

	if (absX > absY)
		return absX * sqrt(1 + pow(absY / absX, 2));
	else
		if (absY == 0)
			return 0;
		else
			return absY * sqrt(1 + pow(absX / absY, 2));
}

double Sign(double A, double b)
{
	return (b > 0 ? fabs(A) : -fabs(A));
}

bool CDiv(double ar, double ai, double br, double bi, double* cr, double* ci)
{
	double s, ars, ais, brs, bis;
	s = fabs(br) + fabs(bi);
	if (s == 0) {
		*cr = 0;
		*ci = 0;
		return false;
	}
	else {
		ars = ar / s;
		ais = ai / s;
		brs = br / s;
		bis = bi / s;
		s = brs * brs + bis * bis;
		*cr = (ars * brs + ais * bis) / s;
		*ci = (ais * brs - ars * bis) / s;
	}
	return true;
}

short sign(double x) {
	return (x > 0) - (x < 0);
}

bool cyno_IsSymmetric(double** theMatrix, long dim1, long dim2, double tol) {
	short i, j;
	double x;
	bool msym;

	msym = false;

	if (dim1 == dim2) {
		for (i = 0; i < dim1; i++)
			for (j = i + 1; j < dim1; j++) {
				x = fabs(theMatrix[i][j] - theMatrix[j][i]);
				if (x > tol)
					return false;
				else if (x > 0 && !msym)
					msym = true;
			}
		if (msym) {
			for (i = 0; i < dim1; i++)
				for (j = i + 1; j < dim1; j++) {
					x = (theMatrix[i][j] + theMatrix[j][i]) / 2;
					theMatrix[i][j] = x;
					theMatrix[j][i] = x;
				}
		}
		return true;
	}
	else
		return false;
}

bool cyno_HouseHolder(long dim, double** aMatrix, double* dVec, double* evec, bool eigenVecs)
{
	short i, j, k, l;
	double sc, hh, h, g, f;

	for (i = dim - 1; i >= 1; i--)
	{
		l = i - 1;
		h = 0;
		sc = 0;
		if (l > 0) {
			for (k = 0; k <= l; k++)
				sc = sc + fabs(aMatrix[i][k]);
			if (sc == 0)
				evec[i] = aMatrix[i][l];
			else {
				for (k = 0; k <= l; k++) {
					aMatrix[i][k] = aMatrix[i][k] / sc;
					h += aMatrix[i][k] * aMatrix[i][k];
				}
				f = aMatrix[i][l];
				if (f >= 0)
					g = -sqrt(h);
				else
					g = sqrt(h);
				evec[i] = sc * g;
				h = h - f * g;
				aMatrix[i][l] = f - g;
				f = 0;
				for (j = 0; j <= l; j++) {
					if (eigenVecs)
						aMatrix[j][i] = aMatrix[i][j] / h;
					g = 0;
					for (k = 0; k <= j; k++)
						g += aMatrix[j][k] * aMatrix[i][k];
					for (k = j + 1; k <= l; k++)
						g += aMatrix[k][j] * aMatrix[i][k];
					evec[j] = g / h;
					f += evec[j] * aMatrix[i][j];
				}
				hh = f / (h + h);
				for (j = 0; j <= l; j++) {
					f = aMatrix[i][j];
					g = evec[j] - hh * f;
					evec[j] = g;
					for (k = 0; k <= j; k++)
						aMatrix[j][k] = aMatrix[j][k] - (f * evec[k] + g * aMatrix[i][k]);
				}
			}
		}
		else
			evec[i] = aMatrix[i][l];
		dVec[i] = h;
	}
	if (eigenVecs)
		dVec[0] = 0;
	evec[0] = 0;
	if (eigenVecs)
		for (i = 0; i < dim; i++) {
			l = i - 1;
			if (dVec[i] != 0)
				for (j = 0; j <= l; j++) {
					g = 0;
					for (k = 0; k <= l; k++)
						g += aMatrix[i][k] * aMatrix[k][j];
					for (k = 0; k <= l; k++)
						aMatrix[k][j] = aMatrix[k][j] - g * aMatrix[k][i];
				}
			dVec[i] = aMatrix[i][i];
			aMatrix[i][i] = 1;
			for (j = 0; j <= l; j++) {
				aMatrix[j][i] = 0;
				aMatrix[i][j] = 0;
			}
		}
	else
		dVec[i] = aMatrix[i][i];
	return true;
}


bool cyno_TriDiagEigenSystem(long numElements, double* dVec, double* evec, double** eigenVecs)
{
	short m, l, iter, i, k;
	double s, r, p, g, f, dd, c, b;
	bool eigvs;
	double tolerance;

	tolerance = 0.0000000001;

	eigvs = true;

	for (i = 1; i < numElements; i++)
		evec[i - 1] = evec[i];
	evec[numElements - 1] = 0;

	for (l = 0; l < numElements; l++) {
		iter = 0;
		do {
			for (m = l; m < numElements - 1; m++) {
				dd = fabs(dVec[m]) + fabs(dVec[m + 1]);
				if (fabs(evec[m]) < tolerance)
					break;
			}
			if (m != l) {
				if (iter == 100)
					return false;
				iter += 1;
				g = (dVec[l + 1] - dVec[l]) / (2 * evec[l]);
				r = Pythag(g, 1);
				g = dVec[m] - dVec[l] + evec[l] / (g + Sign(r, g));
				s = 1;
				c = 1;
				p = 0;
				for (i = m - 1; i >= l; i--) {
					f = s * evec[i];
					b = c * evec[i];
					r = Pythag(f, g);
					evec[i + 1] = r;
					if (r == 0) {
						dVec[i + 1] = dVec[i + 1] - p;
						evec[m] = 0;
						break;
					}
					s = f / r;
					c = g / r;
					g = dVec[i + 1] - p;
					r = (dVec[i] - g) * s + 2 * c * b;
					p = s * r;
					dVec[i + 1] = g + p;
					g = c * r - b;
					if (eigvs)
						for (k = 0; k < numElements; k++) {
							f = eigenVecs[k][i + 1];
							eigenVecs[k][i + 1] = s * eigenVecs[k][i] + c * f;
							eigenVecs[k][i] = c * eigenVecs[k][i] - s * f;
						}
				}
				if (!((r == 0) && (i >= l))) {
					dVec[l] = dVec[l] - p;
					evec[l] = g;
					evec[m] = 0;
				}
			}
		} while (m != l);

	}
	return true;

}



bool UnbalanceVectors(long vMatrixDim1, long vMatrixDim2, double** vMatrix, long transDim,
	double* trans, short bL, short bH, bool isLeft) {
	short i, ii, j, k;
	double s;

	if (vMatrixDim1 == vMatrixDim2 && vMatrixDim2 == transDim) {
		for (i = bL; i <= bH; i++) {
			if (isLeft)
				s = 1 / trans[i];
			else
				s = trans[i];
			for (j = 0; j < vMatrixDim2; j++)
				vMatrix[i][j] = vMatrix[i][j] * s;
		}
		for (ii = 0; ii < transDim; ii++) {
			i = ii;
			if (i < bL || i > bH) {
				if (i < bL)
					i = bL - ii;
				k = trans[i];
				if (k != i)
					for (j = 0; j < vMatrixDim2; j++) {
						s = vMatrix[i][j];
						vMatrix[i][j] = vMatrix[k][j];
						vMatrix[k][j] = s;
					}
			}
		}
		return true;
	}
	else
		return false;
}

bool BalanceMatrix(long theMatrixDim1, long theMatrixDim2, double** theMatrix, long transDim,
	double* trans, short* bL, short* bH) {
	short i, j, rput, cput;
	double tmp, c, r, f, g, s;
	bool isol;

	if (theMatrixDim1 != theMatrixDim2)
		return false;

	rput = theMatrixDim1 - 1;
	cput = 0;

	isol = true;

	for (i = theMatrixDim1 - 1; i >= 0; i--) {
		for (j = 0; j <= rput; j++)
			if (i != j && theMatrix[i][j] != 0) {
				isol = false;
				break;
			}
		if (isol) {
			trans[rput] = i;
			if (rput != i) {
				for (j = 0; j <= theMatrixDim1 - 1; j++) {
					tmp = theMatrix[i][j];
					theMatrix[i][j] = theMatrix[rput][j];
					theMatrix[rput][j] = tmp;
				}
				for (j = 0; j <= theMatrixDim1 - 1; j++) {
					tmp = theMatrix[j][i];
					theMatrix[j][i] = theMatrix[j][rput];
					theMatrix[j][rput] = tmp;
				}
			}
			if (rput == 0) {
				*bL = cput;
				*bH = rput;
				return true;
			}
			rput = rput - 1;
			isol = true;
		}
		else
			isol = true;
	}

	for (i = 0; i <= rput; i++) {
		for (j = 0; j <= rput; j++)
			if (i != j && theMatrix[j][i] != 0) {
				isol = false;
				break;
			}
		if (isol) {
			trans[cput] = i;
			if (cput != i) {
				for (j = 0; j <= theMatrixDim1 - 1; j++) {
					tmp = theMatrix[i][j];
					theMatrix[i][j] = theMatrix[cput][j];
					theMatrix[cput][j] = tmp;
				}
				for (j = 0; j <= theMatrixDim1 - 1; j++) {
					tmp = theMatrix[j][cput];
					theMatrix[j][i] = theMatrix[j][cput];
					theMatrix[j][cput] = tmp;
				}
			}
			if (cput == rput) {
				*bL = cput;
				*bH = rput;
				return true;
			}
			cput++;
			isol = true;
		}
		else
			isol = true;
	}
	for (i = cput; i <= rput; i++)
		trans[i] = 1;
	isol = true;
	while (isol) {
		for (i = cput; i <= rput; i++) {
			c = 0;
			r = 0;
			for (j = cput; j <= rput; j++)
				if (j != i) {
					c += fabs(theMatrix[j][i]);
					r += fabs(theMatrix[i][j]);
				}
			if (c != 0 && r != 0) {
				g = r / RADIX;
				f = 1;
				s = c + r;
				while (c < g) {
					f = f * RADIX;
					c = c * RADIX * RADIX;
				}
				g = r * RADIX;
				while (c >= g) {
					f = f / RADIX;
					c = c / RADIX / RADIX;
				}
				if ((c + r) / f < 0.95 * s) {
					g = 1 / f;
					trans[i] = trans[i] * f;
					isol = true;
					for (j = cput; j <= theMatrixDim1 - 1; j++)
						theMatrix[i][j] = theMatrix[i][j] * g;
					for (j = 0; j <= rput; j++)
						theMatrix[j][i] = theMatrix[j][i] * f;
				}
				else
					isol = false;
			}
			else
				isol = false;
		}
	}
	*bL = cput;
	*bH = rput;
	return true;
}

bool Hessenberg(long aMatrixDim1, long aMatrixDim2, double** aMatrix, long transDim, double* trans,
	long zMatrixDim1, long zMatrixDim2, double** zMatrix, short bL, short bH, bool eigvects) {
	short m, j, i, la, kp;
	short mm, mp, kl, bbl, bbh;
	double x, y;

	if (aMatrixDim1 != aMatrixDim2)
		return false;
	for (m = 0; m < aMatrixDim1; m++)
		trans[m] = 0;
	if (bL == 0)
		bbl = 0;
	else
		bbl = bL;
	if (bH == 0)
		bbh = aMatrixDim1 - 1;
	else
		bbh = bH;
	la = bbh - 1;
	kp = bbl + 1;

	if (la < kp)
		return true;
	else {
		for (m = kp; m <= la; m++) {
			mm = m - 1;
			x = 0;
			i = m;
			for (j = m; j <= bbh; j++)
				if (fabs(aMatrix[j][mm]) > fabs(x)) {
					x = aMatrix[j][mm];
					i = j;
				}
			trans[m] = i;
			if (i != m) {
				for (j = mm; j < aMatrixDim1; j++) {
					y = aMatrix[i][j];
					aMatrix[i][j] = aMatrix[m][j];
					aMatrix[m][j] = y;
				}
				for (j = 0; j <= bbh; j++) {
					y = aMatrix[j][i];
					aMatrix[j][i] = aMatrix[j][m];
					aMatrix[j][m] = y;
				}
			}
			if (x != 0) {
				mp = m + 1;
				for (i = mp; i <= bbh; i++) {
					y = aMatrix[i][mm];
					if (y != 0) {
						y = y / x;
						aMatrix[i][mm] = y;
					}
					for (j = m; j < aMatrixDim1; j++)
						aMatrix[i][j] = aMatrix[i][j] - y * aMatrix[m][j];
					for (j = 0; j <= bbh; j++)
						aMatrix[j][m] = aMatrix[j][m] + y * aMatrix[j][i];
				}
			}
		}
	}
	if (eigvects) {
		for (i = 0; i < aMatrixDim1; i++) {
			for (j = 0; j < aMatrixDim1; j++)
				zMatrix[i][j] = 0;
			zMatrix[i][i] = 1;
		}
		kl = bbh - bbl - 1;
		if (kl >= 1)
			for (mm = bbh - 1; mm >= bbl + 1; mm--) {
				mp = mm + 1;
				for (i = mp; i <= bbh; i++)
					zMatrix[i][mm] = aMatrix[i][mm - 1];
				i = trans[mm];
				if (i != mm) {
					for (j = mm; j <= bbh; j++) {
						zMatrix[mm][j] = zMatrix[i][j];
						zMatrix[i][j] = 0;
					}
					zMatrix[i][mm] = 1;
				}
			}
	}
	// Reset lower triangle
	for (i = bbh; i >= bbl; i--)
		for (j = i - 2; j >= bbl; j--)
			aMatrix[i][j] = 0;
	return true;
}

bool HessEigenSystem(long hMatrixDim1, long hMatrixDim2, double** hMatrix,
	long zMatrixDim1, long zMatrixDim2, double** zMatrix,
	long wrDim, double* wr, long wiDim, double* wi,
	short bL, short bH, bool eigvects) {
	short i, j, k, l, m;
	short en, na, enm2, its, mp2;
	short bbl, bbh;
	double norm, tst1, tst2, x, y, zz, p;
	double q, r, s, ra, sa, t, w;
	double vr, vi, thr, thi;
	bool b;



	if (hMatrixDim1 == hMatrixDim2 && (zMatrixDim1 == zMatrixDim2 || !eigvects)
		&& (hMatrixDim2 == zMatrixDim1 || !eigvects)) {

		if (bL == 0)
			bbl = 0;
		else
			bbl = bL;
		if (bH == 0)
			bbh = hMatrixDim1 - 1;
		else
			bbh = bH;

		norm = 0;
		k = 0;
		for (i = 0; i < hMatrixDim1; i++) {
			for (j = i; j < hMatrixDim1; j++)
				norm = norm + fabs(hMatrix[i][j]);
			if (i < bbl || i > bbh) {
				wr[i] = hMatrix[i][i];
				wi[i] = 0;
			}
		}
		en = bbh;
		t = 0;
		while (en >= bbl) {
			its = 0;
			na = en - 1;
			enm2 = en - 2;
			for (l = en; l >= bbl; l--)
				if (l > bL) {
					s = fabs(hMatrix[l - 1][l - 1]) + fabs(hMatrix[l][l]);
					if (s == 0)
						s = norm;
					tst1 = s;
					tst2 = tst1 + fabs(hMatrix[l][l - 1]);
					if (tst2 == tst1)
						break;
				}
			//			if (l == 0)
			if (l == -1)
				l = bbl;
			x = hMatrix[en][en];

			if (l != en) {
				y = hMatrix[na][na];
				w = hMatrix[en][na] * hMatrix[na][en];

				if (l != na) {
					if (its == 30)
						return false;
					else if (its == 10 || its == 20) {
						t = t + x;
						for (i = bbl; i <= en; i++) {
							hMatrix[i][i] = hMatrix[i][i] - x;

						}
						s = fabs(hMatrix[en][na]) + fabs(hMatrix[na][enm2]);
						x = 0.75 * s;
						y = x;
						w = -0.4375 * s * s;
					}
					its++;
					for (m = en - 2; m >= l; m--) {
						zz = hMatrix[m][m];

						r = x - zz;
						s = y - zz;
						p = (r * s - w) / hMatrix[m + 1][m] + hMatrix[m][m + 1];
						q = hMatrix[m + 1][m + 1] - zz - r - s;

						r = hMatrix[m + 2][m + 1];
						s = fabs(p) + fabs(q) + fabs(r);
						p = p / s;
						q = q / s;

						r = r / s;

						if (m != l) {
							tst1 = fabs(p) * (fabs(hMatrix[m - 1][m - 1]) + fabs(zz) + fabs(hMatrix[m + 1][m + 1]));
							tst2 = tst1 + fabs(hMatrix[m][m - 1]) * (fabs(q) + fabs(r));
							if (tst2 == tst1)
								break;
						}
						else
							break;

					}
					mp2 = m + 2;
					for (i = mp2; i <= en; i++) {
						hMatrix[i][i - 2] = 0;

						if (i != mp2) {
							hMatrix[i][i - 3] = 0;

						}
					}
					for (k = m; k <= na; k++) {
						if (k != m) {
							p = hMatrix[k][k - 1];
							q = hMatrix[k + 1][k - 1];
							r = 0;
							if (k != na)
								r = hMatrix[k + 2][k - 1];
							x = fabs(p) + fabs(q) + fabs(r);
							if (x != 0) {
								p = p / x;
								q = q / x;
								r = r / x;
							}
						}

						s = sqrt(p * p + q * q + r * r) * sign(p);

						if (s != 0) {
							if (k == m) {

								if (l != m) {
									hMatrix[k][k - 1] = -hMatrix[k][k - 1];

								}
							}
							else {
								hMatrix[k][k - 1] = -s * x;

							}
							p = p + s;
							x = p / s;
							y = q / s;
							zz = r / s;
							q = q / p;
							r = r / p;

							if (k != na) {
								for (j = k; j < hMatrixDim1; j++) {
									p = hMatrix[k][j] + q * hMatrix[k + 1][j] + r * hMatrix[k + 2][j];


									hMatrix[k][j] = hMatrix[k][j] - p * x;
									hMatrix[k + 1][j] = hMatrix[k + 1][j] - p * y;
									hMatrix[k + 2][j] = hMatrix[k + 2][j] - p * zz;

								}
								if (en < k + 3)
									j = en;
								else
									j = k + 3;
								for (i = 0; i <= j; i++) {
									p = x * hMatrix[i][k] + y * hMatrix[i][k + 1] + zz * hMatrix[i][k + 2];
									hMatrix[i][k] = hMatrix[i][k] - p;
									hMatrix[i][k + 1] = hMatrix[i][k + 1] - p * q;
									hMatrix[i][k + 2] = hMatrix[i][k + 2] - p * r;

								}
								if (eigvects)
									for (i = bbl; i <= bbh; i++) {
										p = x * zMatrix[i][k] + y * zMatrix[i][k + 1] + zz * zMatrix[i][k + 2];
										zMatrix[i][k] = zMatrix[i][k] - p;
										zMatrix[i][k + 1] = zMatrix[i][k + 1] - p * q;
										zMatrix[i][k + 2] = zMatrix[i][k + 2] - p * r;
									}
							}
							else {
								for (j = k; j <= hMatrixDim1 - 1; j++) {
									p = hMatrix[k][j] + q * hMatrix[k + 1][j];
									hMatrix[k][j] = hMatrix[k][j] - p * x;
									hMatrix[k + 1][j] = hMatrix[k + 1][j] - p * y;

								}
								if (en < k + 3)
									j = en;
								else
									j = k + 3;
								for (i = 0; i <= j; i++) {
									p = x * hMatrix[i][k] + y * hMatrix[i][k + 1];
									hMatrix[i][k] = hMatrix[i][k] - p;
									hMatrix[i][k + 1] = hMatrix[i][k + 1] - p * q;

								}
								if (eigvects)
									for (i = bbl; i <= bbh; i++)
									{
										p = x * zMatrix[i][k] + y * zMatrix[i][k + 1];
										zMatrix[i][k] = zMatrix[i][k] - p;
										zMatrix[i][k + 1] = zMatrix[i][k + 1] - p * q;
									}
							}
						}
					}
				}
				else {//Two roots found
					p = (y - x) / 2;
					q = p * p + w;
					zz = sqrt(fabs(q));
					hMatrix[en][en] = x + t;

					x = hMatrix[en][en];
					hMatrix[na][na] = y + t;

					if (q >= 0) {
						zz = p + fabs(zz) * sign(p);
						wr[na] = x + zz;
						wr[en] = wr[na];
						if (zz != 0)
							wr[en] = x - w / zz;
						wi[na] = 0;
						wi[en] = 0;
						x = hMatrix[en][na];
						s = fabs(x) + fabs(zz);
						p = x / s;
						q = zz / s;
						r = sqrt(p * p + q * q);
						p = p / r;
						q = q / r;
						for (j = na; j < hMatrixDim1; j++) {
							zz = hMatrix[na][j];
							hMatrix[na][j] = q * zz + p * hMatrix[en][j];
							hMatrix[en][j] = q * hMatrix[en][j] - p * zz;

						}
						for (i = 0; i <= en; i++) {
							zz = hMatrix[i][na];
							hMatrix[i][na] = q * zz + p * hMatrix[i][en];
							hMatrix[i][en] = q * hMatrix[i][en] - p * zz;

						}
						if (eigvects)
							for (i = bbl; i <= bbh; i++) {
								zz = zMatrix[i][na];
								zMatrix[i][na] = q * zz + p * zMatrix[i][en];
								zMatrix[i][en] = q * zMatrix[i][en] - p * zz;
							}
					}
					else {
						wr[na] = x + p;
						wr[en] = x + p;
						wi[na] = zz;
						wi[en] = -zz;
					}
					en = enm2;
				}
			}
			else {
				hMatrix[en][en] = x + t;

				wr[en] = hMatrix[en][en];
				wi[en] = 0;
				en = na;
			}
		}
		if (eigvects) {

			if (norm == 0)
				return false;

			for (en = hMatrixDim1 - 1; en >= 0; en--) {
				p = wr[en];
				q = wi[en];
				na = en - 1;

				if (q < 0) {
					m = na;
					if (fabs(hMatrix[en][na]) <= fabs(hMatrix[na][en])) {
						//							CDiv 0#, -hMatrix(na, en), hMatrix(na, na) - p, q, thr, thi
						b = CDiv(0, -hMatrix[na][en], hMatrix[na][na] - p, q, &thr, &thi);
						hMatrix[na][na] = thr;
						hMatrix[na][en] = thi;

					}
					else {
						hMatrix[na][na] = q / hMatrix[en][na];
						hMatrix[na][en] = -(hMatrix[en][en] - p) / hMatrix[en][na];

					}
					hMatrix[en][na] = 0;
					hMatrix[en][en] = 1;

					enm2 = na - 1;
					if (enm2 != -1) //enm2 <> lB -1
						for (i = en - 2; i >= 0; i--) {
							w = hMatrix[i][i] - p;
							ra = 0;
							sa = 0;
							for (j = m; j <= en; j++) {
								ra = ra + hMatrix[i][j] * hMatrix[j][na];
								sa = sa + hMatrix[i][j] * hMatrix[j][en];
							}
							if (wi[i] < 0) {
								zz = w;
								r = ra;
								s = sa;
							}
							else {
								m = i;
								if (wi[i] == 0) {
									//CDiv -ra, -sa, w, q, thr, thi
									b = CDiv(-ra, -sa, w, q, &thr, &thi);
									hMatrix[i][na] = thr;
									hMatrix[i][en] = thi;

								}
								else {
									x = hMatrix[i][i + 1];
									y = hMatrix[i + 1][i];
									vr = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i] - q * q;
									vi = (wr[i] - p) * 2 * q;
									if (vr == 0 && vi == 0) {
										tst1 = norm * (fabs(w) + fabs(q) + fabs(x) + fabs(y) + fabs(zz));
										vr = tst1;
										vr = 0.01 * vr;
										tst2 = tst1 + vr;
										while (tst2 > tst1) {
											vr = 0.01 * vr;
											tst2 = tst1 + vr;
										}
									}
									else {
										//CDiv x * r - zz * ra + q * sa, x * s - zz * sa - q * ra, vr, vi, thr, thi
										b = CDiv(x * r - zz * ra + q * sa, x * s - zz * sa - q * ra, vr, vi, &thr, &thi);
										hMatrix[i][na] = thr;
										hMatrix[i][en] = thi;

									}
									if (fabs(x) <= fabs(zz) + fabs(q)) {
										//CDiv -r - y * hMatrix(i, na), -s - y * hMatrix(i, en), zz, q, thr, thi
										b = CDiv(-r - y * hMatrix[i][na], -s - y * hMatrix[i][en], zz, q, &thr, &thi);
										hMatrix[i + 1][na] = thr;
										hMatrix[i + 1][en] = thi;

									}
									else {
										hMatrix[i + 1][na] = (-ra - w * hMatrix[i][na] + q * hMatrix[i][en]) / x;
										hMatrix[i + 1][en] = (-sa - w * hMatrix[i][en] - q * hMatrix[i][na]) / x;

									}
									if (fabs(hMatrix[i][na]) > fabs(hMatrix[i][en]))
										t = fabs(hMatrix[i][na]);
									else
										t = fabs(hMatrix[i][en]);
									if (t != 0) {
										tst1 = t;
										tst2 = tst1 + 1 / tst1;
										if (tst2 <= tst1) {
											for (j = i; j <= en; j++) {
												hMatrix[j][na] = hMatrix[j][na] / t;
												hMatrix[j][en] = hMatrix[j][en] / t;

											}
										}
									}
								}
							}
						}
				}
				else if (q == 0) {
					m = en;
					hMatrix[en][en] = 1;

					if (na >= 0) {
						for (i = en - 1; i >= 0; i--) {
							w = hMatrix[i][i] - p;
							r = 0;
							for (j = m; j <= en; j++)
								r = r + hMatrix[i][j] * hMatrix[j][en];
							if (wi[i] < 0) {
								zz = w;
								s = r;
							}
							else {
								m = i;
								if (wi[i] == 0) {
									t = w;
									if (t == 0) {
										tst1 = norm;
										t = tst1;
										t = 0.01 * t;
										tst2 = norm + t;
										while (tst2 > tst1) {
											t = 0.01 * t;
											tst2 = norm + t;
										}
									}
									hMatrix[i][en] = -r / t;


								}
								else {
									x = hMatrix[i][i + 1];
									y = hMatrix[i + 1][i];
									q = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i];
									t = (x * s - zz * r) / q;
									hMatrix[i][en] = t;

									if (fabs(x) <= fabs(zz)) {
										hMatrix[i + 1][en] = (-s - y * t) / zz;

									}
									else {
										hMatrix[i + 1][en] = (-r - w * t) / x;

									}
								}
								t = fabs(hMatrix[i][en]);
								if (t != 0) {
									tst1 = t;
									tst2 = tst1 + 1 / tst1;
									if (tst2 <= tst1)
										for (j = i; j <= en; j++) {
											hMatrix[j][en] = hMatrix[j][en] / t;

										}
								}
							}
						}
					}
				}
			}
			for (i = 0; i < hMatrixDim1; i++)
				if (i < bbl || i > bbh)
					for (j = i; j < hMatrixDim1; j++)
						zMatrix[i][j] = hMatrix[i][j];
			for (j = hMatrixDim1 - 1; j >= bbl; j--) {
				if (j < bbh)
					m = j;
				else
					m = bbh;
				for (i = bbl; i <= bbh; i++) {
					zz = 0;
					for (k = bL; k <= m; k++)
						zz = zz + zMatrix[i][k] * hMatrix[k][j];
					zMatrix[i][j] = zz;
				}
			}
		}

		return true;
	}
	else
		return false;
}

//Purpose
//		Eigen System
//Parameters
//		tMatrix		input	square matrix
//		vr			output	1D matrix
//		vi			output	1D matrix
//		vMatrix		output	2D matrix
//Return Values
//		bool
//Remarks
//		tMatrix, vr, vi, vMatrix must be sized properly and all dimensions must have the same size

CYNOMATHUTILITY_API bool _stdcall cyno_EigenSystem(SAFEARRAY** tMatrix, SAFEARRAY** vr, SAFEARRAY** vi, SAFEARRAY** vMatrix,
	bool computeVects, double symTol, bool sortVects) {

	short mtop, mbot;
	bool chk;
	double** ptMatrix, ** pvMatrix;
	double* pArrayElements_vi, * btrans, * htrans;

	long i;

	bool ret;

	long tMatrix_D1Elements, tMatrix_D2Elements;

	tMatrix_D1Elements = (*tMatrix)->rgsabound[0].cElements;
	tMatrix_D2Elements = (*tMatrix)->rgsabound[1].cElements;

	if ((*tMatrix)->cDims == 2 && tMatrix_D1Elements == tMatrix_D2Elements) {

		// Allocate memory for tMatrix and vMatrix assuming vMatrix has same dimension as tMatrix
		ptMatrix = cyno_alloc_matrix(tMatrix_D1Elements, tMatrix_D2Elements);
		pvMatrix = cyno_alloc_matrix(tMatrix_D1Elements, tMatrix_D2Elements);
		btrans = cyno_alloc_double_array(tMatrix_D1Elements);
		htrans = cyno_alloc_double_array(tMatrix_D1Elements);
		if (ptMatrix == NULL || pvMatrix == NULL || btrans == NULL || htrans == NULL) {
			cyno_free_matrix(ptMatrix);
			cyno_free_matrix(pvMatrix);
			free(btrans);
			free(htrans);
			return false;
		}

		ret = cyno_ConvertFromSafeArray2D(tMatrix, tMatrix_D1Elements, tMatrix_D2Elements, ptMatrix);
		ret = cyno_IsSymmetric(ptMatrix, tMatrix_D1Elements, tMatrix_D2Elements, symTol);
		if (ret) {


			if (computeVects) {
				chk = cyno_HouseHolder(tMatrix_D1Elements, ptMatrix, (double*)(*vr)->pvData, (double*)(*vi)->pvData, computeVects);
				if (chk) {

					cyno_MatrixCopy(tMatrix_D1Elements, tMatrix_D2Elements, ptMatrix, pvMatrix);

					chk = cyno_TriDiagEigenSystem(tMatrix_D1Elements, (double*)(*vr)->pvData, (double*)(*vi)->pvData, pvMatrix);

				}
				else
					chk = false;
			}
			else {
				chk = cyno_HouseHolder(tMatrix_D1Elements, ptMatrix, (double*)(*vr)->pvData, (double*)(*vi)->pvData, computeVects);
				if (chk) {
					//ret = cyno_ConvertFromSafeArray2D(vMatrix, tMatrix_D1Elements, tMatrix_D2Elements, pvMatrix);
					chk = cyno_TriDiagEigenSystem(tMatrix_D1Elements, (double*)(*vr)->pvData, (double*)(*vi)->pvData, pvMatrix);
				}
			}

			if (chk) {
				pArrayElements_vi = (double*)(*vi)->pvData;
				for (i = 0; i < (*vi)->rgsabound[0].cElements; i++)
					pArrayElements_vi[i] = 0;
			}
		}
		else {
			if (computeVects) {

				chk = BalanceMatrix(tMatrix_D1Elements, tMatrix_D2Elements, ptMatrix, tMatrix_D1Elements, btrans, &mbot, &mtop);



				if (chk)
					chk = Hessenberg(tMatrix_D1Elements, tMatrix_D2Elements, ptMatrix, tMatrix_D1Elements, htrans, tMatrix_D1Elements,
						tMatrix_D2Elements, pvMatrix, mbot, mtop, computeVects);


				if (chk)
					chk = HessEigenSystem(tMatrix_D1Elements, tMatrix_D2Elements, ptMatrix,
						tMatrix_D1Elements, tMatrix_D2Elements, pvMatrix,
						tMatrix_D1Elements, (double*)(*vr)->pvData, tMatrix_D1Elements, (double*)(*vi)->pvData,
						mbot, mtop, computeVects);


				if (chk)
					chk = UnbalanceVectors(tMatrix_D1Elements, tMatrix_D2Elements, pvMatrix, tMatrix_D1Elements, btrans, mbot, mtop, false);

			}
			else {
				chk = BalanceMatrix(tMatrix_D1Elements, tMatrix_D2Elements, ptMatrix, tMatrix_D1Elements, btrans, &mbot, &mtop);



				if (chk)
					chk = Hessenberg(tMatrix_D1Elements, tMatrix_D2Elements, ptMatrix, tMatrix_D1Elements, htrans, tMatrix_D1Elements,
						tMatrix_D2Elements, pvMatrix, mbot, mtop, computeVects);



				if (chk)
					chk = HessEigenSystem(tMatrix_D1Elements, tMatrix_D2Elements, ptMatrix,
						tMatrix_D1Elements, tMatrix_D2Elements, pvMatrix,
						tMatrix_D1Elements, (double*)(*vr)->pvData, tMatrix_D1Elements, (double*)(*vi)->pvData,
						mbot, mtop, computeVects);

			}
		
		}


		ret = cyno_ConvertToSafeArray2D(vMatrix, tMatrix_D1Elements, tMatrix_D2Elements, pvMatrix);
		ret = cyno_ConvertToSafeArray2D(tMatrix, tMatrix_D1Elements, tMatrix_D2Elements, ptMatrix);


		// De-Allocate memory to prevent memory leak
		cyno_free_matrix(ptMatrix);
		cyno_free_matrix(pvMatrix);
		free(btrans);
		free(htrans);


		return chk;
	}
	else
		return false;
}


//Purpose
//		Simpler massage way - flooring eigenvalues to one value(10^-8)
//Parameters
//		CorrMatrix	input/output	correlation matrix
//Return Values
//		bool
//Remarks
//		The result is returned through CorrMatrix		

CYNOMATHUTILITY_API bool _stdcall cyno_PDCorrMatrix(SAFEARRAY** CorrMatrix) {
	SAFEARRAY* e;
	SAFEARRAY* v;
	SAFEARRAY* vi;
	SAFEARRAY* sa_temp4;
	SAFEARRAYBOUND rgsabound1[2];
	SAFEARRAYBOUND rgsabound2[1];


	short i;
	bool chk, ret;

	long CorrMatrix_D1Elements, CorrMatrix_D2Elements;

	double** peMatrix, ** temp3, ** temp4;
	double* pvMatrix;

	e = NULL;
	v = NULL;
	vi = NULL;
	sa_temp4 = NULL;

	chk = false;
	if ((*CorrMatrix)->cDims == 2) {
		CorrMatrix_D1Elements = (*CorrMatrix)->rgsabound[0].cElements;
		CorrMatrix_D2Elements = (*CorrMatrix)->rgsabound[1].cElements;

		if (CorrMatrix_D1Elements == CorrMatrix_D2Elements) {
			rgsabound1[0].lLbound = 0;
			rgsabound1[0].cElements = CorrMatrix_D1Elements;
			rgsabound1[1].lLbound = 0;
			rgsabound1[1].cElements = CorrMatrix_D1Elements;
			rgsabound2[0].lLbound = 0;
			rgsabound2[0].cElements = CorrMatrix_D1Elements;
			e = SafeArrayCreate(VT_R8, 2, rgsabound1);
			sa_temp4 = SafeArrayCreate(VT_R8, 2, rgsabound1);
			v = SafeArrayCreate(VT_R8, 1, rgsabound2);
			vi = SafeArrayCreate(VT_R8, 1, rgsabound2);
			if (e != NULL && v != NULL && vi != NULL && sa_temp4 != NULL) {


				chk = cyno_EigenSystem(CorrMatrix, &v, &vi, &e, true, 0.0, false);
				if (chk) {

					peMatrix = cyno_alloc_matrix(CorrMatrix_D1Elements, CorrMatrix_D2Elements);
					temp3 = cyno_alloc_matrix(CorrMatrix_D1Elements, CorrMatrix_D2Elements);
					temp4 = cyno_alloc_matrix(CorrMatrix_D1Elements, CorrMatrix_D2Elements);
					pvMatrix = cyno_alloc_double_array(CorrMatrix_D1Elements);
					if (peMatrix == NULL || temp3 == NULL || temp4 == NULL || pvMatrix == NULL) {
						cyno_free_matrix(peMatrix);
						cyno_free_matrix(temp3);
						cyno_free_matrix(temp4);
						free(pvMatrix);
						chk = false;
					} else {

					ret = cyno_ConvertFromSafeArray2D(&e, CorrMatrix_D1Elements, CorrMatrix_D2Elements, peMatrix);
					ret = cyno_ConvertFromSafeArray1D(&v, CorrMatrix_D1Elements, pvMatrix);

					for (i = 0; i < CorrMatrix_D1Elements; i++)
						if (pvMatrix[i] <= 0.00000000001)
							pvMatrix[i] = 0.00000000001;
						else
							//v(i) = Math.Sqr(v(i))
							pvMatrix[i] = pvMatrix[i];
					cyno_MatrixCopy(CorrMatrix_D1Elements, CorrMatrix_D2Elements, peMatrix, temp3);


					chk = cyno_DiagonalMatrixMultiplyMatrixVector(CorrMatrix_D1Elements, CorrMatrix_D1Elements, peMatrix, CorrMatrix_D1Elements, pvMatrix);

					chk = cyno_InverseMatrix(CorrMatrix_D1Elements, CorrMatrix_D1Elements, temp3, temp4);

					ret = cyno_ConvertToSafeArray2D(&sa_temp4, CorrMatrix_D1Elements, CorrMatrix_D1Elements, temp4);
					ret = cyno_ConvertToSafeArray2D(&e, CorrMatrix_D1Elements, CorrMatrix_D1Elements, peMatrix);
					chk = cyno_MatrixMultiply_2(&e, &sa_temp4, CorrMatrix);
					}

					cyno_free_matrix(peMatrix);
					cyno_free_matrix(temp3);
					cyno_free_matrix(temp4);
					free(pvMatrix);

				}
				else
					chk = false;
			}
		}
	}
	if (e != NULL) SafeArrayDestroy(e);
	if (sa_temp4 != NULL) SafeArrayDestroy(sa_temp4);
	if (v != NULL) SafeArrayDestroy(v);
	if (vi != NULL) SafeArrayDestroy(vi);

	return chk;

}

//Purpose
//		Massage the matrix so that we can do cholesky decomposition
//Parameters
//		CorrMatrix	input/output	correlation matrix
//Return Values
//		bool
//Remarks
//		The result is returned through CorrMatrix		

CYNOMATHUTILITY_API bool _stdcall cyno_ApproximateCorrMatrix(SAFEARRAY** CorrMatrix) {
	SAFEARRAY* e;
	SAFEARRAY* v;
	SAFEARRAY* vi;
	SAFEARRAY* sa_temp3;
	SAFEARRAY* sa_temp4;
	SAFEARRAYBOUND rgsabound1[2];
	SAFEARRAYBOUND rgsabound2[1];



	short i, j;
	bool chk, ret;
	double s;

	long CorrMatrix_D1Elements, CorrMatrix_D2Elements;

	double** peMatrix, * temp2, ** temp3;
	double* pvMatrix;

	e = NULL;
	v = NULL;
	vi = NULL;
	sa_temp3 = NULL;
	sa_temp4 = NULL;

	chk = false;
	if ((*CorrMatrix)->cDims == 2) {
		CorrMatrix_D1Elements = (*CorrMatrix)->rgsabound[0].cElements;
		CorrMatrix_D2Elements = (*CorrMatrix)->rgsabound[1].cElements;

		if (CorrMatrix_D1Elements == CorrMatrix_D2Elements) {
			rgsabound1[0].lLbound = 0;
			rgsabound1[0].cElements = CorrMatrix_D1Elements;
			rgsabound1[1].lLbound = 0;
			rgsabound1[1].cElements = CorrMatrix_D1Elements;
			rgsabound2[0].lLbound = 0;
			rgsabound2[0].cElements = CorrMatrix_D1Elements;
			e = SafeArrayCreate(VT_R8, 2, rgsabound1);
			sa_temp3 = SafeArrayCreate(VT_R8, 2, rgsabound1);
			sa_temp4 = SafeArrayCreate(VT_R8, 2, rgsabound1);
			v = SafeArrayCreate(VT_R8, 1, rgsabound2);
			vi = SafeArrayCreate(VT_R8, 1, rgsabound2);
			if (e != NULL && v != NULL && vi != NULL && sa_temp4 != NULL && sa_temp3 != NULL) {


				chk = cyno_EigenSystem(CorrMatrix, &v, &vi, &e, true, 0.0, false);
				if (chk) {

					peMatrix = cyno_alloc_matrix(CorrMatrix_D1Elements, CorrMatrix_D2Elements);
					temp3 = cyno_alloc_matrix(CorrMatrix_D1Elements, CorrMatrix_D2Elements);
					pvMatrix = cyno_alloc_double_array(CorrMatrix_D1Elements);
					temp2 = cyno_alloc_double_array(CorrMatrix_D1Elements);
					if (peMatrix == NULL || temp3 == NULL || pvMatrix == NULL || temp2 == NULL) {
						cyno_free_matrix(peMatrix);
						cyno_free_matrix(temp3);
						free(pvMatrix);
						free(temp2);
						chk = false;
					} else {

					ret = cyno_ConvertFromSafeArray2D(&e, CorrMatrix_D1Elements, CorrMatrix_D2Elements, peMatrix);
					ret = cyno_ConvertFromSafeArray1D(&v, CorrMatrix_D1Elements, pvMatrix);

					for (i = 0; i < CorrMatrix_D1Elements; i++)
						if (pvMatrix[i] < 0)
							pvMatrix[i] = 0;
						else
							pvMatrix[i] = sqrt(pvMatrix[i]);

					for (i = 0; i < CorrMatrix_D1Elements; i++) {
						s = 0;
						for (j = 0; j < CorrMatrix_D1Elements; j++)
							s = s + peMatrix[i][j] * peMatrix[i][j] * pvMatrix[j];
						temp2[i] = sqrt(1 / s);
					}

					chk = cyno_DiagonalMatrixMultiplyMatrixVector(CorrMatrix_D1Elements, CorrMatrix_D1Elements, peMatrix, CorrMatrix_D1Elements, pvMatrix);


					chk = cyno_DiagonalMatrixMultiply(CorrMatrix_D1Elements, temp2, CorrMatrix_D1Elements, CorrMatrix_D1Elements, peMatrix);


					cyno_MatrixCopy(CorrMatrix_D1Elements, CorrMatrix_D2Elements, peMatrix, temp3);



					ret = cyno_ConvertToSafeArray2D(&e, CorrMatrix_D1Elements, CorrMatrix_D1Elements, peMatrix);
					chk = cyno_MatrixTranspose(&e, &sa_temp4);
					ret = cyno_ConvertFromSafeArray2D(&sa_temp4, CorrMatrix_D1Elements, CorrMatrix_D2Elements, peMatrix);



					ret = cyno_ConvertToSafeArray2D(&sa_temp3, CorrMatrix_D1Elements, CorrMatrix_D1Elements, temp3);
					chk = cyno_MatrixMultiply_2(&sa_temp3, &sa_temp4, CorrMatrix);
					}

					cyno_free_matrix(peMatrix);
					cyno_free_matrix(temp3);
					free(pvMatrix);
					free(temp2);

				}
				else
					chk = false;
			}
		}
	}
	if (e != NULL) SafeArrayDestroy(e);
	if (sa_temp3 != NULL) SafeArrayDestroy(sa_temp3);
	if (sa_temp4 != NULL) SafeArrayDestroy(sa_temp4);
	if (v != NULL) SafeArrayDestroy(v);
	if (vi != NULL) SafeArrayDestroy(vi);

	return chk;

}

//Purpose
//		Solving a set of linear equations (simple case, only for non-singular matrix)
//Parameters
//		xInputJacob		input	square matrix
//		xInputArray		input	1D matrix
//		xOutputArray	output	1D matrix
//Return Values
//		bool
//Remarks
//		All dimensions must be in same size. xInputJacob, xInputArray, xOutPutArray must be sized properly before call this function.		

CYNOMATHUTILITY_API bool _stdcall cyno_LinearSystemSolver(SAFEARRAY** xInputJacob, SAFEARRAY** xInputArray, SAFEARRAY** xOutputArray)
{
	double** temp, ** inverse;
	long D1Elements, D2Elements;
	bool chk;
	long i;



	chk = false;
	if ((*xInputJacob)->cDims == 2) {
		D1Elements = (*xInputJacob)->rgsabound[0].cElements;
		D2Elements = (*xInputJacob)->rgsabound[1].cElements;
		if (D1Elements == D2Elements) {
			// allocate memory
			temp = cyno_alloc_matrix(D1Elements, D2Elements);
			inverse = cyno_alloc_matrix(D1Elements, D2Elements);
			if (temp == NULL || inverse == NULL) {
				cyno_free_matrix(temp);
				cyno_free_matrix(inverse);
				return false;
			}

			chk = cyno_ConvertFromSafeArray2D(xInputJacob, D1Elements, D2Elements, temp);
			chk = cyno_InverseMatrix(D1Elements, D2Elements, temp, inverse);
			chk = cyno_ConvertToSafeArray2D(xInputJacob, D1Elements, D1Elements, inverse);
			chk = cyno_MatrixMultiply_2(xInputJacob, xInputArray, xOutputArray);

			cyno_free_matrix(temp);
			cyno_free_matrix(inverse);
		}
		else
			chk = false;
	}

	return chk;
}

CYNOMATHUTILITY_API bool __stdcall cyno_MatrixScalePlain(double a, const double* input, long rows, long cols, double* output)
{
	long i;
	long size;

	if (input == NULL || output == NULL || rows <= 0 || cols <= 0)
		return false;

	size = rows * cols;
	for (i = 0; i < size; i++)
		output[i] = a * input[i];

	return true;
}

CYNOMATHUTILITY_API bool __stdcall cyno_MatrixMultiplyPlain(const double* a, long aRows, long aCols, const double* b, long bRows, long bCols, double* output)
{
	double** aMatrix;
	double** bMatrix;
	double** outMatrix;
	bool ret;

	if (a == NULL || b == NULL || output == NULL || aRows <= 0 || aCols <= 0 || bRows <= 0 || bCols <= 0)
		return false;

	aMatrix = cyno_alloc_matrix(aRows, aCols);
	bMatrix = cyno_alloc_matrix(bRows, bCols);
	outMatrix = cyno_alloc_matrix(aRows, bCols);
	if (aMatrix == NULL || bMatrix == NULL || outMatrix == NULL) {
		cyno_free_matrix(aMatrix);
		cyno_free_matrix(bMatrix);
		cyno_free_matrix(outMatrix);
		return false;
	}

	cyno_copy_flat_to_matrix(a, aRows, aCols, aMatrix);
	cyno_copy_flat_to_matrix(b, bRows, bCols, bMatrix);
	ret = cyno_MatrixMulti2D2D(aRows, aCols, aMatrix, bRows, bCols, bMatrix, outMatrix);
	if (ret)
		cyno_copy_matrix_to_flat(outMatrix, aRows, bCols, output);

	cyno_free_matrix(aMatrix);
	cyno_free_matrix(bMatrix);
	cyno_free_matrix(outMatrix);
	return ret;
}

CYNOMATHUTILITY_API bool __stdcall cyno_VectorDotPlain(const double* a, const double* b, long size, double* output)
{
	long i;
	double value;

	if (a == NULL || b == NULL || output == NULL || size <= 0)
		return false;

	value = 0.0;
	for (i = 0; i < size; i++)
		value += a[i] * b[i];
	*output = value;
	return true;
}

CYNOMATHUTILITY_API bool __stdcall cyno_InverseMatrixPlain(const double* input, long dim, double* output)
{
	double** inMatrix;
	double** outMatrix;
	bool ret;

	if (input == NULL || output == NULL || dim <= 0)
		return false;

	inMatrix = cyno_alloc_matrix(dim, dim);
	outMatrix = cyno_alloc_matrix(dim, dim);
	if (inMatrix == NULL || outMatrix == NULL) {
		cyno_free_matrix(inMatrix);
		cyno_free_matrix(outMatrix);
		return false;
	}

	cyno_copy_flat_to_matrix(input, dim, dim, inMatrix);
	ret = cyno_InverseMatrix(dim, dim, inMatrix, outMatrix);
	if (ret)
		cyno_copy_matrix_to_flat(outMatrix, dim, dim, output);

	cyno_free_matrix(inMatrix);
	cyno_free_matrix(outMatrix);
	return ret;
}

CYNOMATHUTILITY_API bool __stdcall cyno_MatrixTransposePlain(const double* input, long rows, long cols, double* output)
{
	double** inMatrix;
	double** outMatrix;

	if (input == NULL || output == NULL || rows <= 0 || cols <= 0)
		return false;

	inMatrix = cyno_alloc_matrix(rows, cols);
	outMatrix = cyno_alloc_matrix(cols, rows);
	if (inMatrix == NULL || outMatrix == NULL) {
		cyno_free_matrix(inMatrix);
		cyno_free_matrix(outMatrix);
		return false;
	}

	cyno_copy_flat_to_matrix(input, rows, cols, inMatrix);
	cyno_MatrixTrans(rows, cols, inMatrix, outMatrix);
	cyno_copy_matrix_to_flat(outMatrix, cols, rows, output);

	cyno_free_matrix(inMatrix);
	cyno_free_matrix(outMatrix);
	return true;
}

CYNOMATHUTILITY_API bool __stdcall cyno_MatrixPowerPlain(double* inOutMatrix, long dim, short power)
{
	double** inMatrix;
	double** tmpMatrix;
	double** outMatrix;
	long i, j;
	short rPow;
	bool ret;

	if (inOutMatrix == NULL || dim <= 0 || power < 1)
		return false;

	inMatrix = cyno_alloc_matrix(dim, dim);
	tmpMatrix = cyno_alloc_matrix(dim, dim);
	outMatrix = cyno_alloc_matrix(dim, dim);
	if (inMatrix == NULL || tmpMatrix == NULL || outMatrix == NULL) {
		cyno_free_matrix(inMatrix);
		cyno_free_matrix(tmpMatrix);
		cyno_free_matrix(outMatrix);
		return false;
	}

	cyno_copy_flat_to_matrix(inOutMatrix, dim, dim, inMatrix);
	cyno_MatrixCopy(dim, dim, inMatrix, tmpMatrix);
	rPow = 1;
	ret = true;

	while (ret && rPow < power) {
		ret = cyno_MatrixMulti2D2D(dim, dim, tmpMatrix, dim, dim, inMatrix, outMatrix);
		rPow++;
		cyno_MatrixCopy(dim, dim, outMatrix, tmpMatrix);
	}

	if (ret) {
		if (power == 1)
			cyno_copy_matrix_to_flat(inMatrix, dim, dim, inOutMatrix);
		else
			cyno_copy_matrix_to_flat(outMatrix, dim, dim, inOutMatrix);
	}

	cyno_free_matrix(inMatrix);
	cyno_free_matrix(tmpMatrix);
	cyno_free_matrix(outMatrix);
	return ret;
}

CYNOMATHUTILITY_API bool __stdcall cyno_CholeskyDecompPlain(double* inOutMatrix, long dim)
{
	double** matrix;
	double** tmp;
	short i, j, k;
	double sum;

	if (inOutMatrix == NULL || dim <= 0)
		return false;

	matrix = cyno_alloc_matrix(dim, dim);
	tmp = cyno_alloc_matrix(dim, dim);
	if (matrix == NULL || tmp == NULL) {
		cyno_free_matrix(matrix);
		cyno_free_matrix(tmp);
		return false;
	}

	cyno_copy_flat_to_matrix(inOutMatrix, dim, dim, matrix);
	for (i = 0; i < dim; i++)
		for (j = i; j < dim; j++) {
			sum = matrix[i][j];
			for (k = i - 1; k >= 0; k--)
				sum = sum - matrix[i][k] * matrix[j][k];
			if (i == j) {
				if (sum < 0) {
					cyno_free_matrix(matrix);
					cyno_free_matrix(tmp);
					return false;
				}
				tmp[i][i] = sqrt(sum);
			}
			else {
				tmp[j][i] = sum / tmp[i][i];
				matrix[j][i] = tmp[j][i];
			}
		}
	for (i = 0; i < dim; i++)
		for (j = i + 1; j < dim; j++)
			tmp[i][j] = 0;

	cyno_copy_matrix_to_flat(tmp, dim, dim, inOutMatrix);
	cyno_free_matrix(matrix);
	cyno_free_matrix(tmp);
	return true;
}

CYNOMATHUTILITY_API bool __stdcall cyno_EigenSystemPlain(const double* input, long dim, double* vr, double* vi, double* eigenVectors, bool computeVects, double symTol, bool sortVects)
{
	short mtop, mbot;
	bool chk;
	bool ret;
	double** tMatrix;
	double** vMatrix;
	double* btrans;
	double* htrans;
	long i;

	(void)sortVects;

	if (input == NULL || vr == NULL || vi == NULL || eigenVectors == NULL || dim <= 0)
		return false;

	tMatrix = cyno_alloc_matrix(dim, dim);
	vMatrix = cyno_alloc_matrix(dim, dim);
	btrans = cyno_alloc_double_array(dim);
	htrans = cyno_alloc_double_array(dim);
	if (tMatrix == NULL || vMatrix == NULL || btrans == NULL || htrans == NULL) {
		cyno_free_matrix(tMatrix);
		cyno_free_matrix(vMatrix);
		free(btrans);
		free(htrans);
		return false;
	}

	cyno_copy_flat_to_matrix(input, dim, dim, tMatrix);
	ret = cyno_IsSymmetric(tMatrix, dim, dim, symTol);
	if (ret) {
		chk = cyno_HouseHolder(dim, tMatrix, vr, vi, computeVects);
		if (chk) {
			if (computeVects)
				cyno_MatrixCopy(dim, dim, tMatrix, vMatrix);
			chk = cyno_TriDiagEigenSystem(dim, vr, vi, vMatrix);
		}
		if (chk) {
			for (i = 0; i < dim; i++)
				vi[i] = 0.0;
		}
	} else {
		chk = BalanceMatrix(dim, dim, tMatrix, dim, btrans, &mbot, &mtop);
		if (chk)
			chk = Hessenberg(dim, dim, tMatrix, dim, htrans, dim, dim, vMatrix, mbot, mtop, computeVects);
		if (chk)
			chk = HessEigenSystem(dim, dim, tMatrix, dim, dim, vMatrix, dim, vr, dim, vi, mbot, mtop, computeVects);
		if (chk && computeVects)
			chk = UnbalanceVectors(dim, dim, vMatrix, dim, btrans, mbot, mtop, false);
	}

	if (chk)
		cyno_copy_matrix_to_flat(vMatrix, dim, dim, eigenVectors);

	cyno_free_matrix(tMatrix);
	cyno_free_matrix(vMatrix);
	free(btrans);
	free(htrans);
	return chk;
}

CYNOMATHUTILITY_API bool __stdcall cyno_PDCorrMatrixPlain(double* inOutMatrix, long dim)
{
	double* eigenValues;
	double* imagValues;
	double* eigenVectors;
	double** peMatrix;
	double** temp3;
	double** temp4;
	bool chk;
	long i;

	if (inOutMatrix == NULL || dim <= 0)
		return false;

	eigenValues = cyno_alloc_double_array(dim);
	imagValues = cyno_alloc_double_array(dim);
	eigenVectors = cyno_alloc_double_array(dim * dim);
	peMatrix = cyno_alloc_matrix(dim, dim);
	temp3 = cyno_alloc_matrix(dim, dim);
	temp4 = cyno_alloc_matrix(dim, dim);
	if (eigenValues == NULL || imagValues == NULL || eigenVectors == NULL || peMatrix == NULL || temp3 == NULL || temp4 == NULL) {
		free(eigenValues);
		free(imagValues);
		free(eigenVectors);
		cyno_free_matrix(peMatrix);
		cyno_free_matrix(temp3);
		cyno_free_matrix(temp4);
		return false;
	}

	chk = cyno_EigenSystemPlain(inOutMatrix, dim, eigenValues, imagValues, eigenVectors, true, 0.0, false);
	if (chk) {
		cyno_copy_flat_to_matrix(eigenVectors, dim, dim, peMatrix);
		for (i = 0; i < dim; i++)
			if (eigenValues[i] <= 1e-11)
				eigenValues[i] = 1e-11;
		cyno_MatrixCopy(dim, dim, peMatrix, temp3);
		chk = cyno_DiagonalMatrixMultiplyMatrixVector(dim, dim, peMatrix, dim, eigenValues);
		if (chk)
			chk = cyno_InverseMatrix(dim, dim, temp3, temp4);
		if (chk)
			chk = cyno_MatrixMulti2D2D(dim, dim, peMatrix, dim, dim, temp4, temp3);
		if (chk)
			cyno_copy_matrix_to_flat(temp3, dim, dim, inOutMatrix);
	}

	free(eigenValues);
	free(imagValues);
	free(eigenVectors);
	cyno_free_matrix(peMatrix);
	cyno_free_matrix(temp3);
	cyno_free_matrix(temp4);
	return chk;
}

CYNOMATHUTILITY_API bool __stdcall cyno_ApproximateCorrMatrixPlain(double* inOutMatrix, long dim)
{
	double* eigenValues;
	double* imagValues;
	double* eigenVectors;
	double* temp2;
	double** peMatrix;
	double** temp3;
	double** temp4;
	bool chk;
	long i, j;
	double s;

	if (inOutMatrix == NULL || dim <= 0)
		return false;

	eigenValues = cyno_alloc_double_array(dim);
	imagValues = cyno_alloc_double_array(dim);
	eigenVectors = cyno_alloc_double_array(dim * dim);
	temp2 = cyno_alloc_double_array(dim);
	peMatrix = cyno_alloc_matrix(dim, dim);
	temp3 = cyno_alloc_matrix(dim, dim);
	temp4 = cyno_alloc_matrix(dim, dim);
	if (eigenValues == NULL || imagValues == NULL || eigenVectors == NULL || temp2 == NULL || peMatrix == NULL || temp3 == NULL || temp4 == NULL) {
		free(eigenValues);
		free(imagValues);
		free(eigenVectors);
		free(temp2);
		cyno_free_matrix(peMatrix);
		cyno_free_matrix(temp3);
		cyno_free_matrix(temp4);
		return false;
	}

	chk = cyno_EigenSystemPlain(inOutMatrix, dim, eigenValues, imagValues, eigenVectors, true, 0.0, false);
	if (chk) {
		cyno_copy_flat_to_matrix(eigenVectors, dim, dim, peMatrix);
		for (i = 0; i < dim; i++)
			eigenValues[i] = (eigenValues[i] < 0.0) ? 0.0 : sqrt(eigenValues[i]);

		for (i = 0; i < dim; i++) {
			s = 0.0;
			for (j = 0; j < dim; j++)
				s += peMatrix[i][j] * peMatrix[i][j] * eigenValues[j];
			temp2[i] = sqrt(1.0 / s);
		}

		chk = cyno_DiagonalMatrixMultiplyMatrixVector(dim, dim, peMatrix, dim, eigenValues);
		if (chk)
			chk = cyno_DiagonalMatrixMultiply(dim, temp2, dim, dim, peMatrix);
		if (chk) {
			cyno_MatrixCopy(dim, dim, peMatrix, temp3);
			cyno_MatrixTrans(dim, dim, peMatrix, temp4);
			chk = cyno_MatrixMulti2D2D(dim, dim, temp3, dim, dim, temp4, peMatrix);
		}
		if (chk)
			cyno_copy_matrix_to_flat(peMatrix, dim, dim, inOutMatrix);
	}

	free(eigenValues);
	free(imagValues);
	free(eigenVectors);
	free(temp2);
	cyno_free_matrix(peMatrix);
	cyno_free_matrix(temp3);
	cyno_free_matrix(temp4);
	return chk;
}

CYNOMATHUTILITY_API bool __stdcall cyno_LinearSystemSolverPlain(const double* jacobian, const double* rhs, long dim, double* output)
{
	double** temp;
	double** inverse;
	double* result;
	bool chk;
	long i;

	if (jacobian == NULL || rhs == NULL || output == NULL || dim <= 0)
		return false;

	temp = cyno_alloc_matrix(dim, dim);
	inverse = cyno_alloc_matrix(dim, dim);
	result = cyno_alloc_double_array(dim);
	if (temp == NULL || inverse == NULL || result == NULL) {
		cyno_free_matrix(temp);
		cyno_free_matrix(inverse);
		free(result);
		return false;
	}

	cyno_copy_flat_to_matrix(jacobian, dim, dim, temp);
	chk = cyno_InverseMatrix(dim, dim, temp, inverse);
	if (chk)
		chk = cyno_MatrixMulti2D1D(dim, dim, inverse, dim, (double*)rhs, result);
	if (chk)
		for (i = 0; i < dim; i++)
			output[i] = result[i];

	cyno_free_matrix(temp);
	cyno_free_matrix(inverse);
	free(result);
	return chk;
}
