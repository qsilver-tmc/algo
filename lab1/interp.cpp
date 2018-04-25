// ����� ������� ��� ������������ 
#include "interp.h"
#include <stdlib.h>
#include <stdio.h>

/* ������������ ���������
 * x - �����, f(x) ������� ���� �����
 * arrX - ������ �����, �� ������� ���������� ������������. (���������������)
 * koefsTable - ������� ������� ������������ �������������
 * return �������� ������� � ����e (x)
*/
double interp_splain(double const x, std::vector<double> arrX, std::vector<koefs_t> koefsTable)
{
	int pos = find_pos(x, arrX.data(), arrX.size()) + 1;
	double dx = x - arrX[pos-1];
    return koefsTable[pos].a + koefsTable[pos].b*dx + koefsTable[pos].c*dx*dx + koefsTable[pos].d*dx*dx*dx;
}

/* ���������� ������� �������������
 * arrX, arrValues - ������ �����, �� ������� ���������� ������������. (���������������)
 * koefsTable - ������� ��� ������ ��������� (������ �������� �������)
*/
void fillKoefsTable (std::vector<double> arrX, std::vector<double> arrValues, 
					   std::vector<koefs_t> &koefsTable)
{
    double A, B, D, F;
    koefsTable[1].ksi = koefsTable[1].eta = 0;
    for (size_t i = 2; i < arrX.size(); i++)
    {
        A = arrX[i-1] - arrX[i-2];
        B = -2 * (arrX[i] - arrX[i-2]);
        D = arrX[i] - arrX[i-1];
        F = (-3) * ((arrValues[i] - arrValues[i-1]) / D - (arrValues[i-1] - arrValues[i-2]) / A);

        koefsTable[i+1].ksi = D / (B - A * koefsTable[i].ksi);
        koefsTable[i+1].eta = (F + A * koefsTable[i].eta)/(B - A * koefsTable[i].ksi);
    }

    koefsTable[arrX.size()-1].c = 0;
	double dx;
    for (int i = arrX.size() - 1; i > 0; i--)
    {
        dx = arrX[i] - arrX[i-1];
        koefsTable[i].c = koefsTable[i+1].ksi * koefsTable[i+1].c + koefsTable[i+1].eta;
		koefsTable[i].a = arrValues[i-1];
        koefsTable[i].d = (koefsTable[i+1].c - koefsTable[i].c)/(3*dx);
        koefsTable[i].b = (arrValues[i] - arrValues[i-1])/dx - dx*(koefsTable[i+1].c + 2*koefsTable[i].c)/3;
    }
}

/* ����������� ������������ ��������� �������
 * (x, y) - �����, f(x,y) ������� ���� �����
 * degreeX - ������� �������� �� �, degreeY - ������� �������� �� Y
 * arrX, arrY, tableValues - ������� �����, �� ������� ���������� ������������. (���������������)
 * return �������� ������� � ����e (x,y)
*/
double interp_newton_3d(double const x, double const y, int degreeX, int degreeY, 
					 std::vector<double> arrX, std::vector<double> arrY, 
					 std::vector<std::vector<double> > tableValues)
{
	double funcXY = 0; // ������� ��������
	int nearestX, nearestY; // ����� �����, � ������� ���������� ������ ��������� ������ degree + 1
	nearestX = search_nearest(degreeX + 1, x, arrX.data(), arrX.size());
	nearestY = search_nearest(degreeY + 1, y, arrY.data(), arrY.size());
	
	std::vector<double> arrFirstInterp (degreeX + 1);
	for(int i = 0; i < degreeX + 1; i++)
	{
		arrFirstInterp[i] = 
		calc_div_dif_matr(&(arrY.data()[nearestY]), &(tableValues.data()[nearestX + i].data()[nearestY]), degreeY + 1, y);
	}
	funcXY = calc_div_dif_matr(&(arrX.data()[nearestX]), arrFirstInterp.data(), degreeX + 1, x);
	return funcXY;
}

/* ������������ ��������� �������
 * x - �����, Y ������� ���� �����
 * degree - ������� ��������
 * arrX, arrY - ������� �����, �� ������� ���������� ������������. (���������������)
 * size ������� ��������
 * return �������� ������� � ����e X
*/
double interp_newton(double const x, int degree, double *arrX, double *arrY, int const size)
{
	double y = 0; // ������� ��������
	int nearest; // ����� �����, � ������� ���������� ������ ��������� � ������� ������ degree + 1

	nearest = search_nearest(degree + 1, x, arrX, size);

#ifdef _DEBUG
		int i;
		printf("\n");
		for (i = nearest; i < nearest + degree + 1; i++)
			printf("\n[debug] %f", arrX[i]);
		printf("\n\n");
#endif

	y = calc_div_dif_matr(&(arrX[nearest]), &(arrY[nearest]), degree+1, x); 

	return y;
}

/* ����� ��������� ���������
 * count ���������� ������ ��������� �����
 * destPoint �����, � ������� ������ ���������
 * dataBase ������ � ������� (������ ���� ������������)
 * size ����� �������
 * return startEnum ����� �������� �������, � �������� ���������� ������������ ���������
*/
int search_nearest(int count, double point, double *arr, int size)
{
	int startEnum = 0;
	if (point < arr[0]) // �����������. �����
		startEnum =  0;
	else if (point > arr[size - 1]) // ����������� ������
		startEnum = size - count;
	else 
	{
		// ������� ���� ����� � arr, ����� ������� ���� X
		startEnum = find_pos(point, arr, size);

		if (startEnum < count/2)
			startEnum =  0;
		else if (startEnum > size - ((count % 2 == 0) ? count/2 : (count/2 + 1)) - 1)
			return size - count;
		else
			startEnum -= count/2 ;
	}
	return startEnum;
}

/* ������������� �������� ����������� ���������
 * col1 - ��������� �� ������ ������� ����������� ��������
 * col2 ��������� �� ������
 * size ����� ��������
 * dest �����, ��� ������� ���� ������ ��������
 * return ��������� ����������
*/
double calc_div_dif_matr(double *col1, double *col2, int size, double dest)
{
	int i, j; // ��������
	double **divDif = 0; // ������� ����������� ���������
	double iter; // ������������ ���������� ��� ��������
	double result = 0;
	
	//[cols][rows]
	divDif = (double **)realloc(divDif, sizeof(double *)*(size+1) + sizeof(double)*summ(1, size-1));
	assert(divDif != 0);
	divDif[0] = col1;
	divDif[1] = col2;
	for (i = 2; i <= size+1; i++)
		divDif[i] = (double*)((char*)divDif + (size+1)*sizeof(double*) + summ(size+2-i, size-1)*sizeof(double));

	for (i = 2; i <= size; i++) // �������� �� ��������, ������� � ����������� �������� �� ���� ��-���
		for (j = 0; j < size - i + 1; j++) // �������� �� �������, ��������� j-� ��-� i-���� �������
			divDif[i][j] = (divDif[i-1][j] - divDif[i-1][j+1]) / (divDif[0][j] - divDif[0][j + (i - 1)]);

	// ������� �������� ��������
	iter = 1;
	for (i = 0; i < size; i++)
	{
		result += iter * divDif[i + 1][0];
		iter *= (dest - divDif[0][i]);
	}

	free(divDif);
	return result;
}

// ���������� ����� ����� �� A �� B (������������)
int summ(int a, int b)
{
	int sum = 0;
	for (int i = a; i < b + 1; i++)
		sum += i;
	return sum;
}

// ����� ������� - ������� ��������������� �������, �����, � �������� ���������� - ��������
int find_pos(double valueForSearch, double *arr, int size)
{
	int left = 0, right = 0;
	int &leftRef = left, &rightRef = right;
	log_find(valueForSearch, 0, arr, size, leftRef, rightRef);

	return bin_find(valueForSearch, leftRef, rightRef, arr);
}

// ��������������� ����� � arr ������� � startIndex. ��� ������ ������������� left � right
int log_find(double valueForSearch, int startIndex, double *arr, int size, int &left, int &right)
{
	int offset = 1;
	for(;;)
	{
		offset = offset << 1;
		if(startIndex + offset >= size)  // ����� �� �������
		{
			right = size - 1;
			left = startIndex + (offset >> 1) - 1;
			break;
		}
		else if(arr[startIndex + offset - 1] > valueForSearch) // ����� ������ �������
		{
			right = startIndex + offset;
			left = startIndex + (offset >> 1) - 1; 
			break;
		}
	}

	return 0;
}

// �������� ����� � arr ����� first � last
int bin_find(double valueForSearch, int left, int right, double *arr)
{
	int mid = 0;
	while(left < right)
	{
		mid = (left + right) / 2; // �� �� ����������� ������������ �� �����
		if(valueForSearch < arr[mid] || mid == left)
			right = mid;
		else
			left = mid;
	}
	return mid; 
}