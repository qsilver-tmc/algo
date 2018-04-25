// набор функций дл€ интерпол€ции 
#include "interp.h"
#include <stdlib.h>
#include <stdio.h>

/* »нтерпол€ци€ сплайнами
 * x - точка, f(x) которой надо найти
 * arrX - массив точек, по которым проводитс€ »нтерпол€ци€. (отсортированный)
 * koefsTable - таблица заранее подсчитанных коэффициентов
 * return значение функции в точкe (x)
*/
double interp_splain(double const x, std::vector<double> arrX, std::vector<koefs_t> koefsTable)
{
	int pos = find_pos(x, arrX.data(), arrX.size()) + 1;
	double dx = x - arrX[pos-1];
    return koefsTable[pos].a + koefsTable[pos].b*dx + koefsTable[pos].c*dx*dx + koefsTable[pos].d*dx*dx*dx;
}

/* «аполнение таблицы коэффициентов
 * arrX, arrValues - массив точек, по которым проводитс€ »нтерпол€ци€. (отсортированный)
 * koefsTable - таблица дл€ записи результат (ѕјћя“№ ¬џƒ≈Ћ≈Ќј «ј–јЌ≈≈)
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

/* ћЌќ√ќћ≈–Ќјя »нтерпол€ци€ полиномом Ќьютона
 * (x, y) - точка, f(x,y) которой надо найти
 * degreeX - степень полинома по ’, degreeY - степень полинома по Y
 * arrX, arrY, tableValues - массивы точек, по которым проводитс€ »нтерпол€ци€. (отсортированный)
 * return значение функции в точкe (x,y)
*/
double interp_newton_3d(double const x, double const y, int degreeX, int degreeY, 
					 std::vector<double> arrX, std::vector<double> arrY, 
					 std::vector<std::vector<double> > tableValues)
{
	double funcXY = 0; // искомое значение
	int nearestX, nearestY; // номер точки, с которой начинаетс€ список ближайших длиной degree + 1
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

/* »нтерпол€ци€ полиномом Ќьютона
 * x - точка, Y которой надо найти
 * degree - степень полинома
 * arrX, arrY - массивы точек, по которым проводитс€ »нтерпол€ци€. (отсортированный)
 * size размеры массивов
 * return значение функции в точкe X
*/
double interp_newton(double const x, int degree, double *arrX, double *arrY, int const size)
{
	double y = 0; // искомое значение
	int nearest; // номер точки, с которой начинаетс€ список ближайших к искомой длиной degree + 1

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

/* ѕоиск ближайших элементов
 * count количество нужных ближайших точек
 * destPoint точка, к которой ищутс€ ближайшие
 * dataBase массив с точками (должен быть отсортирован)
 * size длина массива
 * return startEnum номер элемента массива, с которого начинаетс€ перечисление ближайших
*/
int search_nearest(int count, double point, double *arr, int size)
{
	int startEnum = 0;
	if (point < arr[0]) // экстрапол€ц. снизу
		startEnum =  0;
	else if (point > arr[size - 1]) // экстрапол€ц сверху
		startEnum = size - count;
	else 
	{
		// сначала ищем точку в arr, после которой идет X
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

/* »нтерполирует матрицей разделенных разностей
 * col1 - указатель на первый столбец разделенной разности
 * col2 указатель на второй
 * size длина массивов
 * dest точка, дл€ которой надо искать значение
 * return результат вычислений
*/
double calc_div_dif_matr(double *col1, double *col2, int size, double dest)
{
	int i, j; // итерации
	double **divDif = 0; // матрица разделенных разностей
	double iter; // итерационна€ переменна€ дл€ подсчета
	double result = 0;
	
	//[cols][rows]
	divDif = (double **)realloc(divDif, sizeof(double *)*(size+1) + sizeof(double)*summ(1, size-1));
	assert(divDif != 0);
	divDif[0] = col1;
	divDif[1] = col2;
	for (i = 2; i <= size+1; i++)
		divDif[i] = (double*)((char*)divDif + (size+1)*sizeof(double*) + summ(size+2-i, size-1)*sizeof(double));

	for (i = 2; i <= size; i++) // итераци€ по столбцам, начина€ с разделенной разности от двух эл-тов
		for (j = 0; j < size - i + 1; j++) // итераци€ по строкам, вычисл€ем j-й эл-т i-того столбца
			divDif[i][j] = (divDif[i-1][j] - divDif[i-1][j+1]) / (divDif[0][j] - divDif[0][j + (i - 1)]);

	// —„»“ј≈ћ »“ќ√ќ¬ќ≈ «Ќј„≈Ќ»≈
	iter = 1;
	for (i = 0; i < size; i++)
	{
		result += iter * divDif[i + 1][0];
		iter *= (dest - divDif[0][i]);
	}

	free(divDif);
	return result;
}

// вычисление суммы чисел от A до B (включительно)
int summ(int a, int b)
{
	int sum = 0;
	for (int i = a; i < b + 1; i++)
		sum += i;
	return sum;
}

// поиск позиции - сначала логарифмическим поиском, потом, в найденом промежутке - бинарным
int find_pos(double valueForSearch, double *arr, int size)
{
	int left = 0, right = 0;
	int &leftRef = left, &rightRef = right;
	log_find(valueForSearch, 0, arr, size, leftRef, rightRef);

	return bin_find(valueForSearch, leftRef, rightRef, arr);
}

// логарифмический поиск в arr начина€ с startIndex. ѕри выходе устанавливает left и right
int log_find(double valueForSearch, int startIndex, double *arr, int size, int &left, int &right)
{
	int offset = 1;
	for(;;)
	{
		offset = offset << 1;
		if(startIndex + offset >= size)  // выход за пределы
		{
			right = size - 1;
			left = startIndex + (offset >> 1) - 1;
			break;
		}
		else if(arr[startIndex + offset - 1] > valueForSearch) // нашли нужный отрезок
		{
			right = startIndex + offset;
			left = startIndex + (offset >> 1) - 1; 
			break;
		}
	}

	return 0;
}

// бинарный поиск в arr между first и last
int bin_find(double valueForSearch, int left, int right, double *arr)
{
	int mid = 0;
	while(left < right)
	{
		mid = (left + right) / 2; // из за ограничений переполнени€ не боюсь
		if(valueForSearch < arr[mid] || mid == left)
			right = mid;
		else
			left = mid;
	}
	return mid; 
}