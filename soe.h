#pragma once

#ifndef SOE_H
#define SOE_H

#include <iostream>
#include "grid.h"

struct Soe {
	double** H;
	double* P;

	Soe(Grid grid) {

		H = new double* [grid.nN];
		P = new double[grid.nN];

		for (int z = 0; z < grid.nN; z++) {
			H[z] = new double[grid.nN];

			P[z] = 0.0;
			for (int j = 0; j < grid.nN; j++) {
				H[z][j] = 0.0;
			}
		}

		for (int i = 0; i < grid.nE; i++) {
			for (int j = 0; j < 4; j++) {
				for (int k = 0; k < 4; k++) {
					H[grid.elements[i].ID[j]-1][grid.elements[i].ID[k]-1] += grid.elements[i].H[j][k];
				}

				P[grid.elements[i].ID[j]-1] += grid.elements[i].P[j];
			}
		}

		for (int i = 0; i < grid.nN; i++) {
			for (int j = 0; j < grid.nN; j++) {
				cout << H[j][j] << " ";
			}
			cout << endl;
		}

		for (int i = 0; i < grid.nN; i++) {
			cout << P[i] << " ";
		}
		cout << endl;

	}

	void ukladRownanGauss(Grid grid, double** M) {
        int size = grid.nN;

        double m = 0; int k = 0;
        for (int i = 0; i < size; i++) {

            for (int z = i; z < size - 1; z++) {

                if (M[i][i] != 0) {
                    m = M[z + 1][i] / M[i][i];

                    for (int j = i; j < size + 1; j++) {

                        M[z + 1][j] -= M[i][j] * m;

                    }
                }
                else {
                    cout << "0 na przekatnej " << endl;
                    break;
                }
            }

        }
        //wypisanie macierzy
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size + 1; j++) {
                cout << M[i][j] << "\t";
            }
            cout << "\n";
        }
        cout << endl;
        double x = 0.0, suma = 0.0;
        double* tabx = new double[size - 1];
        tabx[size - 1] = M[size - 1][size] / M[size - 1][size - 1];
        for (int i = size - 2; i >= 0; i--) {
            suma = 0;
            for (int j = i + 1; j < size; j++) {
                suma += M[i][j] * tabx[j];
            }
            tabx[i] = (M[i][size] - suma) / M[i][i];

        }

        for (int i = 0; i < size; i++) {
            cout << "Rozwiazanie ukladu rownan: x" << i << "=" << tabx[i] << endl;
        }
	}

};


#endif
