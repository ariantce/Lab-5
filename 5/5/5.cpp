#include <iostream>
#include <cmath>

double F1(double x)
{
    return pow((x + pow(x, 3)), 0.5);
}

double F2(double x)
{
    return pow(pow(x, 3) - 1, -0.5);
}

double f(double x, double y)
{
    return (pow(x, 2)) / (1 + pow(y, 2));
}

double trapezoid_rule(double (*f) (double), double a, double b, int nSeg = 1, double e = 1e-8)
{
    double dx, sum, res = 0, resPrev;
    do {
        resPrev = res;
        dx = 1 * (b - a) / nSeg;      //1.0
        sum = 0.5 * (f(a) + f(b));
        for (int i = 1; i < nSeg; i++) {
            sum += f(a + i * dx);
        }
        res = sum * dx;
        nSeg *= 2;
    } while (fabs(res - resPrev) > 3 * e);
    return res;
}

double Simpson_rule(double (*f) (double), double a, double b, int nSeg = 2, double e = 1e-8)
{
    if (nSeg % 2 != 0) {
        std::cout << "THE NUMBER OF NODES MUST BE EVEN.\n";
        exit(0);
    }
    double dx, sum, res = 0, resPrev;
    do {
        resPrev = res;
        dx = (b - a) / nSeg;
        sum = f(a) + f(b);
        for (int i = 1; i < nSeg; i += 2) {
            sum += f(a + i * dx) * 4;
        }
        for (int i = 2; i < nSeg; i += 2) {
            sum += f(a + i * dx) * 2;
        }
        res = dx / 3 * sum;
        nSeg *= 2;
    } while (fabs(res - resPrev) > 15 * e);
    return res;
}

double Simpson_rule(double (*f) (double, double), double a, double b, double c, double d, int m = 2, int n = 2, double e = 1e-8) {
    if (n % 2 != 0 || m % 2 != 0) {
        std::cout << "THE NUMBER OF NODES MUST BE EVEN.\n";
        exit(0);
    }
    double dx = (b - a) / (2 * n);
    double dy = (d - c) / (2 * m);
    double sum = 0;
    double res = 0;
    double resPrev;

    do {
        resPrev = res;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                sum += f(a + (2 * i * dx), c + (2 * j * dy));
                sum += 4 * f(a + ((2 * i + 1) * dx), c + (2 * j * dy));
                sum += f(a + ((2 * i + 2) * dx), c + (2 * j * dy));
                sum += 4 * f(a + (2 * i * dx), c + ((2 * j + 1) * dy));
                sum += 16 * f(a + ((2 * i + 1) * dx), c + ((2 * j + 1) * dy));
                sum += 4 * f(a + ((2 * i + 2) * dx), c + ((2 * j + 1) * dy));
                sum += f(a + (2 * i * dx), c + ((2 * j + 2) * dy));
                sum += 4 * f(a + ((2 * i + 1) * dx), c + ((2 * j + 2) * dy));
                sum += f(a + ((2 * i + 2) * dx), c + ((2 * j + 2) * dy));
            }
        }
        res = sum * (dx * dy / 9);
        sum = 0;
    } while (fabs(res - resPrev) > 15 * e);
    return res;
}



int main() {
    double res = Simpson_rule(f, 0, 4, 1, 2, 50, 50);
    std::cout << res << '\n';
    return 0;
}