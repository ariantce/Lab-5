#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace std;

const double E1 = 1e-4;
const double E2 = 1e-5;

double F1(double x)
{
	return (pow((x + x * x * x), 0.5));
}

double F2(const double& x, const double& y)
{
	return (x * x + 2 * y);
}

double Simpson(double a, double b, double e, double F1(double))
{
	int n = 2;
	double h = (b - a) / n;
	double I1 = F1(a) + F1(b);
	double I2;
	double x_i = a;
	double sum_Nech = 0;
	double sum_Ch = 0;

	do {

		I2 = I1;
		n *= 2;
		h = (b - a) / n;
		x_i = a + h;

		for (int i = 0; i < n; i++) {
			if ((i % 2) == 0)
			{
				sum_Ch += F1(x_i);
			}
			else {
				sum_Nech += F1(x_i);
			}
			x_i += h;
		}

		I1 = (h / 3) * (F1(a) + 4 * sum_Nech + 2 * sum_Ch + F1(b));
	} while (fabs(I1 - I2) > 15 * e);

	return I1;
}

double Simpson_Kub(const double& a, const double& b, const double& c, const double& d, double func2(const double&, const double&), const double& N, const double& M)
{

	double h_x = (b - a) / (2 * N);
	double h_y = (d - c) / (2 * M);

	double sum = 0;
	double I = 0;
	double X_i = a;
	double Y_i = c;

	vector<double>X;
	vector<double>Y;

	do {
		X.push_back(X_i);
		X_i += h_x;
	} while (X_i <= b);

	do {
		Y.push_back(Y_i);
		Y_i += h_y;
	} while (Y_i <= d);

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			sum += func2(X[2 * i], Y[2 * j]) + 4 * func2(X[2 * i + 1], Y[2 * j]) +
				func2(X[2 * i + 2], Y[2 * j]) + 4 * func2(X[2 * i], Y[2 * j + 1]) +
				16 * func2(X[2 * i + 1], Y[2 * j + 1]) + 4 * func2(X[2 * i + 2], Y[2 * j + 1])
				+ func2(X[2 * i], Y[2 * j + 2]) + 4 * func2(X[2 * i + 1], Y[2 * j + 2])
				+ func2(X[2 * i + 2], Y[2 * j + 2]);
		}
	}
	I += sum;
	I *= ((h_x * h_y) / 9);

	return I;
}

double Trapezoid(double a, double b, const double E, double func1(double)) {

	int n = 2;
	double h = (b - a) / n;
	double sum = 0;
	double I1 = func1(a) + func1(b);
	double I2, X_n;

	do {
		I2 = I1;
		n *= 2;
		h = (b - a) / n;
		X_n = a;
		sum = 0;
		for (int i = 0; i < n - 1; i++) 
		{
			X_n += h;
			sum += func1(X_n);
		}
		I1 = (h / 2) * (func1(a) + 2 * sum + func1(b));
	} while (fabs(I1 - I2) >= 3 * E);

	return I1;
}

int main() {+
	double a1 = 0.6;
	double b1 = 1.724;
	double a2 = 0;
	double b2 = 2.0;
	double c = 0;
	double d = 1.0;

	cout << "Trapezoid's method: " << endl;
	cout << setprecision(8) << Trapezoid(a1, b1, E1, F1) << endl;
	cout << "Simson's method: " << endl;
	cout << setprecision(8) << Simpson(a1, b1, E2, F1) << endl;

	int n, m;
	cout << "enter N and M:" << endl;
	cin >> n;
	cin >> m;

	cout << "Simpson's cubature method: " << endl;
	cout << setprecision(8) << Simpson_Kub(a2, b2, c, d, F2, n, m) << endl;

	return 0;
}