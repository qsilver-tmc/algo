#pragma once

#include <iostream> // ��� ������ ���������� ����
#include <math.h>
#include <assert.h>
#include <vector>

#define _DEBUG

typedef struct koefs{
    double a;
    double b;
    double c;
    double d;
    double ksi;
    double eta;
} koefs_t;

/* ������������ ���������
 * x - �����, f(x) ������� ���� �����
 * arrX, arrY - ������� �����, �� ������� ���������� ������������. (���������������)
 * return �������� ������� � ����e (x)
*/
double interp_splain(double const x, std::vector<double> arrX, std::vector<koefs_t> koefsTable);

/* ���������� ������� �������������
 * arrX, arrValues - ������ �����, �� ������� ���������� ������������. (���������������)
 * koefsTable - ������� ��� ������ ��������� (������ �������� �������)
*/
void fillKoefsTable (std::vector<double> arrX, std::vector<double> arrValues, 
					   std::vector<koefs_t> &koefsTable);

/* ����������� ������������ ��������� �������
 * (x, y) - �����, f(x,y) ������� ���� �����
 * degreeX - ������� �������� �� �, degreeY - ������� �������� �� Y
 * arrX, arrY, tableValues - ������� �����, �� ������� ���������� ������������. (���������������)
 * return �������� ������� � ����e (x,y)
*/
double interp_newton_3d(double const x, double const y, int degreeX, int degreeY, 
					 std::vector<double> arrX, std::vector<double> arrY, 
					 std::vector<std::vector<double> > tableValues);

/* ������������ ��������� �������
 * x - �����, Y ������� ���� �����
 * degree - ������� ��������
 * arrX, arrY - ������� �����, �� ������� ���������� ������������. (���������������)
 * size ������� ��������
 * return �������� ������� � ����e X
*/
double interp_newton(double const x, int degree, double *arrX, double *arrY, int const size);

/* ����� ��������� ���������
 * count ���������� ������ ��������� �����
 * destPoint �����, � ������� ������ ���������
 * dataBase ������ � ������� (������ ���� ������������)
 * size ����� �������
 * return startEnum ����� �������� �������, � �������� ���������� ������������ ���������
*/
int search_nearest(int count, double destPoint, double *dataBase, int size);

/* ������������� �������� ����������� ���������
 * col1 - ��������� �� ������ ������� ����������� ��������
 * col2 ��������� �� ������
 * size ����� ��������
 * dest �����, ��� ������� ���� ������ ��������
 * return ��������� ����������
*/
double calc_div_dif_matr(double *col1, double *col2, int size, double dest);

// ���������� ����� ����� �� A �� B (������������)
int summ(int a, int b);

// ����� ������� - ������� ��������������� �������, �����, � �������� ���������� - ��������
int find_pos(double valueForSearch, double *arr, int size);

// ��������������� ����� � arr ������� � startIndex. ��� ������ ������������� left � right
int log_find(double valueForSearch, int startIndex, double *arr, int size, int &left, int &right);

// �������� ����� � arr ����� first � last
int bin_find(double valueForSearch, int left, int right, double *arr);