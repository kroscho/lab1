#include <windows.h>
#include "iostream"
#include "cmath"
#include <consoleapi2.h>
#include <iomanip>


using namespace std;

const int var = 3;
const int n = 6;
double b = 2.0;
double h = 1 / 5.;

//исходна€ функци€
double F(double x)
{
	if (var == 3)
		return pow(5, x) - 3;
	else
		return pow(3, x) - 2 * x + 5;
}

//производна€ исходной функции
double DF(double x)
{
	if (var == 3)
		return log(5) * pow(5, x);
	else
		return log(3) * pow(3, x) - 2;
}

//получение значений x и f(x)
void get_X_Y_Arr(double* xArr, double* xArr1, double* yArr)
{
	double a = 1.0;
	for (int i = 0; i < n; i++)
	{
		xArr[i] = 1.0 + i * h;
	}
	for (int i = 0; i < n - 1; i++)
	{
		xArr1[i] = 1.0 + (i + 0.5) * h;
	}
	for (int i = 0; i < n; i++)
	{
		yArr[i] = F(xArr[i]);
	}
}

//получение таблицы разделенных разностей
void FindDividedDifferenceMatrix(double** DivDiff, double* xArr, double* yArr)
{
	for (int i = 0; i < n; i++)
	{
		DivDiff[i][0] = xArr[i];
		DivDiff[i][1] = yArr[i];
	}
	for (int j = 2; j <= n; j++)
	{
		for (int i = 0; i < n - j + 1; i++)
			DivDiff[i][j] = (DivDiff[i + 1][j - 1] - DivDiff[i][j - 1]) / (xArr[i + j - 1] - xArr[i]);
	}
}

double fi0(double tau)
{
	return (1 + 2 * tau) * pow(1 - tau, 2);
}

double fi1(double tau)
{
	return tau * pow(1 - tau, 2);
}

// погрешность интерпол€ции сплайнами
double FindErrorForSpline(double M4, double M5)
{
	return (M4 / 384 + M5 * h / 240) * pow(h, 4);
}

//метод кубических сплайнов дефекта 1
void SplineInterpolation(double* xArr1, double* xArr2)
{
	cout << "\n»нтерпол€ци€ кубическим сплайном\n";

	double* m = new double[n];
	double df0 = DF(1), dfn = DF(2);

	m[0] = df0;
	m[n - 1] = dfn;

	// находим параметры m методом прогонки
	double* alpha = new double[n];
	double* beta = new double[n];

	for (int j = 1; j < n - 1; j++)
		m[j] = (F(xArr1[j + 1]) - F(xArr1[j - 1])) / (2 * h);

	/*alpha[0] = -1/4;
	beta[0] = 3/(4*h*h)*(F(xArr1[2])-2*F(xArr1[1])+F(xArr1[0]));

	for (int j = 2; j < n; j++)
	{
		alpha[j - 1] = -1 / (4 + alpha[j - 2]);
		beta[j - 1] = (3/(h*h)*(F(xArr1[0]) - 2 * F(xArr1[j-1]) + F(xArr1[j-2])-beta[j-2]))/(4+alpha[j-2]);
	}

	for (int j = n - 1; j >= 1; j--)
		m[j - 1] = alpha[j - 1] * m[j] + beta[j - 1];*/
	
	//тут верно
	/*alpha[1] = 0;
	beta[1] = df0;

	for (int j = 1; j < n - 1; j++)
	{
		alpha[j + 1] = -1 / (4 + alpha[j]);
		beta[j + 1] = (3 * (F(xArr1[j + 1]) - F(xArr1[j - 1])) / h - beta[j]) / (4 + alpha[j]);
	}

	for (int j = n - 2; j >= 0; j--)
		m[j] = alpha[j + 1] * m[j + 1] + beta[j + 1];*/
	//до сюда верно

	// наибольшее значение 4-й и 5-й производной достигаетс€ в точке х = 2
	double M4, M5;
	M4 = (var == 3) ? pow(log(5), 4) * 25 : pow(log(3), 4) * 9;
	M5 = (var == 3) ? pow(log(5), 4) * 25 : pow(log(3), 5) * 9;

	cout << "M5 = " << M5 << "\nM4 = " << M4 << endl << endl;

	// оценка погрешности аппроксимации 1-й производной
	cout << " x[i] | " << " df/dx(x[i]) " << "  m[i]   " << " delta " << " оценка " << endl;
	for (int i = 0; i < n; i++)
	{
		double DF_xi = DF(xArr1[i]);
		cout << setprecision(1) << xArr1[i] << "   | " << setprecision(6) << DF_xi << " | " << m[i] << " | " <<
			abs(DF_xi - m[i]) << " | " << M5 / 60 * pow(h, 4) << endl;
	}

	// оценка погрешности аппроксимации функции
	cout << endl << " x  | " << "   f(x)   " << "   S31(f;x) " << " Abs(f(x)-S31(f;x)) " << " оценка " << endl;
	for (int i = 0; i < n - 1; i++)
	{
		double tau = 0.5; // (xArr2[i] - xArr1[i]) / h;
		double S = fi0(tau) * F(xArr1[i]) + fi0(1 - tau) * F(xArr1[i + 1]) +
			h * (fi1(tau) * m[i] - fi1(1 - tau) * m[i + 1]);

		double F_ = F(xArr2[i]);

		cout << setprecision(1) << xArr2[i] << "  | " << setprecision(6) << F_ << " | " <<
			S << " | "  << abs(F_ - S) << " | " <<  FindErrorForSpline(M4, M5) << endl;
	}

	delete[] m;
	delete[] alpha;
	delete[] beta;
}

//печать таблицы разделенных разностей
void PrintDividedDifferenceMatrix(double **DivDiff)
{
	cout << "“аблица разделенных разностей:\n";
	for (int i = 0; i < n; i++)
	{		
		for (int j = 0; j <= n - i; j++)		
			cout << fixed << setprecision(6) << setw(12) << DivDiff[i][j];
		cout << endl;
	}
}

// погрешность интерпол€ции по формуле Ќьютона
double FindErrorForNewton(double x, double* xArr, double M6)
{
	double w = 1;
	for (int i = 0; i < n; i++)
		w *= x - xArr[i];
	int fact = 1;
	for (int i = 2; i < n + 1; i++)
		fact *= i;

	return M6 * abs(w) / fact;
}

//метод Ќьютона
void Newton(double *xArr, double *xArr1, double *yArr)
{
	cout << "»нтерпол€ционна€ формула Ќьютона\n";

	double** DivDiff = new double* [n];

	for (int i = 0; i < n; i++)
	{
		DivDiff[i] = new double[n + 1];
	}

	get_X_Y_Arr(xArr, xArr1, yArr);

	FindDividedDifferenceMatrix(DivDiff, xArr, yArr);

	PrintDividedDifferenceMatrix(DivDiff);

	double M6 = (var == 3) ? 25 * pow(log(5), 6) : 9 * pow(log(3), 6);
	cout << "\nM6 = " << M6 << endl;

	cout << "\nx    |    f(x)  |    Pn(x) |   Delta  |   ќценка |\n";

	for (int k = 0; k < n - 1; k++)
	{
		double f = DivDiff[0][1];
		double w = 1;
		for (int i = 1; i < n; i++)
		{
			w *= xArr1[k] - xArr[i - 1];
			f += DivDiff[0][i + 1] * w;
		}
		double F_ = F(xArr1[k]);
		cout << setprecision(1) << xArr1[k] << "  | " << setprecision(6) << F_ << " | " <<
			f << " | " << abs(f - F_) << " | " << FindErrorForNewton(xArr1[k], xArr, M6) << " | \n";
	}
	for (int i = 0; i < n; i++)
		delete[] DivDiff[i];
	delete[] DivDiff;
}

int main()
{
	//задаем кодировку дл€ вывода символов на экран
	SetConsoleCP(1251);
	//задаем кодировку дл€ ввода символов с клавиатуры в консоль
	SetConsoleOutputCP(1251);

	double* xArr = new double[n];
	double* xArr1 = new double[n];
	double* yArr = new double[n];

	get_X_Y_Arr(xArr, xArr1, yArr);

	Newton(xArr, xArr1, yArr);

	SplineInterpolation(xArr, xArr1);



	cout << endl;
	delete[] xArr;
	delete[] xArr1;
	delete[] yArr;
	system("pause");
}