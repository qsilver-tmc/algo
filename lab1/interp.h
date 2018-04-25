#pragma once

#include <iostream> // дл€ вывода отладочной инфы
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

/* »нтерпол€ци€ сплайнами
 * x - точка, f(x) которой надо найти
 * arrX, arrY - массивы точек, по которым проводитс€ »нтерпол€ци€. (отсортированный)
 * return значение функции в точкe (x)
*/
double interp_splain(double const x, std::vector<double> arrX, std::vector<koefs_t> koefsTable);

/* «аполнение таблицы коэффициентов
 * arrX, arrValues - массив точек, по которым проводитс€ »нтерпол€ци€. (отсортированный)
 * koefsTable - таблица дл€ записи результат (ѕјћя“№ ¬џƒ≈Ћ≈Ќј «ј–јЌ≈≈)
*/
void fillKoefsTable (std::vector<double> arrX, std::vector<double> arrValues, 
					   std::vector<koefs_t> &koefsTable);

/* ћЌќ√ќћ≈–Ќјя »нтерпол€ци€ полиномом Ќьютона
 * (x, y) - точка, f(x,y) которой надо найти
 * degreeX - степень полинома по ’, degreeY - степень полинома по Y
 * arrX, arrY, tableValues - массивы точек, по которым проводитс€ »нтерпол€ци€. (отсортированный)
 * return значение функции в точкe (x,y)
*/
double interp_newton_3d(double const x, double const y, int degreeX, int degreeY, 
					 std::vector<double> arrX, std::vector<double> arrY, 
					 std::vector<std::vector<double> > tableValues);

/* »нтерпол€ци€ полиномом Ќьютона
 * x - точка, Y которой надо найти
 * degree - степень полинома
 * arrX, arrY - массивы точек, по которым проводитс€ »нтерпол€ци€. (отсортированный)
 * size размеры массивов
 * return значение функции в точкe X
*/
double interp_newton(double const x, int degree, double *arrX, double *arrY, int const size);

/* ѕоиск ближайших элементов
 * count количество нужных ближайших точек
 * destPoint точка, к которой ищутс€ ближайшие
 * dataBase массив с точками (должен быть отсортирован)
 * size длина массива
 * return startEnum номер элемента массива, с которого начинаетс€ перечисление ближайших
*/
int search_nearest(int count, double destPoint, double *dataBase, int size);

/* »нтерполирует матрицей разделенных разностей
 * col1 - указатель на первый столбец разделенной разности
 * col2 указатель на второй
 * size длина массивов
 * dest точка, дл€ которой надо искать значение
 * return результат вычислений
*/
double calc_div_dif_matr(double *col1, double *col2, int size, double dest);

// вычисление суммы чисел от A до B (включительно)
int summ(int a, int b);

// поиск позиции - сначала логарифмическим поиском, потом, в найденом промежутке - бинарным
int find_pos(double valueForSearch, double *arr, int size);

// логарифмический поиск в arr начина€ с startIndex. ѕри выходе устанавливает left и right
int log_find(double valueForSearch, int startIndex, double *arr, int size, int &left, int &right);

// бинарный поиск в arr между first и last
int bin_find(double valueForSearch, int left, int right, double *arr);