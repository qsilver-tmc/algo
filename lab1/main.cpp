#include "interp.h"
#include <math.h>
#include <iomanip>

#define PI 3.14159265

double func(double x) {
	return sin ((30*x)*PI/180);
}

using namespace std;

int main() {
	//------ GIVEN -------------
	double x 	= 1.3; 	// x   |
	int n 		= 2; 	// n   |
	int size 	= 7;	// size|
	//--------------------------

	double* arrX = new double[size];
	double* arrY = new double[size];

	for (int i = 0; i < size; i++)
	{
		arrX[i] = -3 + i;
		arrY[i] = func(arrX[i]);
	}

	cout << endl << "interp: " << setprecision(5) << interp_newton(x, n, arrX, arrY, size) << endl;
	cout << "real:   " << func(x) << endl;

	//delete[] arrX, arrY;
	return 0;
}