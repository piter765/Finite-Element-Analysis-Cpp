#pragma once
#ifndef METODA_GAUSSA_H
#define METODA_GAUSSA_H

#include <iostream>
#include <string>
#include <cmath>

using namespace std;

struct Pw {
	double* p;
	double* w;

	Pw(int size) {
		p = new double[size];
		w = new double[size];
	}
};

double f1(double x) {
	return 5*x*x + 3*x + 6;
}

double f2(double x, double y) {
	return 5 * x * x * y * y + 3 * x * y + 6;
}


Pw metodaGaussa(int przestrzen, int points) {
	Pw pw = Pw(points);

	if (points == 2) {
		pw.p[0] = -sqrt(1.0 / 3.0);
		pw.p[1] = sqrt(1.0 / 3.0);

		pw.w[0] = 1.0;
		pw.w[1] = 1.0;
	}
	else if (points == 3) {
		pw.p[0] = -sqrt(3.0 / 5.0);
		pw.p[1] = 0.0;
		pw.p[2] = sqrt(3.0 / 5.0);

		pw.w[0] = 5.0 / 9.0;
		pw.w[1] = 8.0 / 9.0;
		pw.w[2] = 5.0 / 9.0;
	}
	else if (points == 4) {
		pw.p[0] = -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));
		pw.p[1] = -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));
		pw.p[2] = sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));
		pw.p[3] = sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));

		pw.w[0] = (18.0 + sqrt(30.0)) / 36.0;
		pw.w[1] = (18.0 - sqrt(30.0)) / 36.0;
		pw.w[2] = (18.0 + sqrt(30.0)) / 36.0;
		pw.w[3] = (18.0 - sqrt(30.0)) / 36.0;
	}

	double suma = 0;
	for (int i = 0; i < points; i++) {
		if (przestrzen == 2) {
			suma += pw.w[i] * f1(pw.p[i]);
		}
		else if (przestrzen == 3) {
			for (int j = 0; j < points; j++) {
				suma += pw.w[i] * pw.w[j] * f2(pw.p[i], pw.p[j]);
			}
		}
	}

	return pw;

}


#endif 