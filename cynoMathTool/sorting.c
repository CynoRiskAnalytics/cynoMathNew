#include "pch.h"
#include "framework.h"
#include "math.h"
#include <stdlib.h>

//Purpose
//		Heap sorting
//Parameters
//		*x		The first element of a vector to be sorted
//		*index	The first element of a vector which stores the sort result
//		*iMin	Lower bound of the vector to be sorted 
//		*IMax	Upper bound of the vector to be sorted
//Return Values
//		integer
//Remarks
//		Vectors x and index must be sized properly in caller
//		Lower bound of x and index must be 0 and
//		it should be called like cyno_sort(index(0), x(0), ...)
//Change history
//		Jan 31, 2020 - created from an older version by Jesse

CYNOMATHUTILITY_API int __stdcall cyno_Sort(int* index1, double* x, int n)
{
	int i, ir, j, l, irra;
	int m;
	int* index;

	if (index1 == NULL || x == NULL || n < 0) {
		return -1;
	}

	if (n < 2) {
		if (n == 1) {
			index1[0] = 0;
		}
		return 0;
	}

	m = n + 1;
	index = (int*)malloc((size_t)m * sizeof(int));
	if (index == NULL) {
		return -1;
	}
	
	for (i = 1; i <= n; i++)
	{
		index[i] = i - 1;
	}

	l = (int)(n / 2.0 + 1 + 0.5);
	ir = n;

	while (1) {
		if (l > 1) {
//			l = l - 1;
			irra = index[--l];
		}
		else {
			irra = index[ir];
			index[ir] = index[1];
//			ir = ir - 1;
			if (--ir == 1) {
				index[1] = irra;
				break;
			}
		}

		i = l;
		j = l + l;

		while (j <= ir) {
			if (j < ir)
			{
				if (x[index[j]] < x[index[j + 1]])
					j = j + 1;
			}
			if (x[irra] < x[index[j]]) {
				index[i] = index[j];
				i = j;
				j = j + j;
			}
			else
				break;
		}
		index[i] = irra;

	}

	for (i = 0; i < n; i++) {
		index1[i] = index[i + 1];
	}

	free(index);

	return 0;
}

