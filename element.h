#pragma once
#ifndef ELEMENT_H
#define ELEMENT_H

#include "grid.h"
#include "elementUniwersalny.h"

struct Node {
	double x;
	double y;
	double t;
	int BC;

	Node(double x1 = 0.0, double y1 = 0.0, double t1 = 0.0, int BC = 0) : x(x1), y(y1), t(t1), BC(0) {};
};


struct Element {
	int ID[4];
	double H[4][4];
	double HBC[4][4];
	int liczbaPunktowIntegracji = 2;
	double P[4];

	Element() : liczbaPunktowIntegracji(2) {
		for (int i = 0; i < 4; i++) {
			ID[i] = 0;
			for (int j = 0; j < 4; j++) {
				H[i][j] = 0.0;
				HBC[i][j] = 0.0;
				P[i] = 0.0;
			}
		}
	}

	void createElement(int liczbaPunktow, Node* nodes, double conductivity, double alfa, double tot) {
		liczbaPunktowIntegracji = liczbaPunktow;

		ElementUniwersalny element = ElementUniwersalny(liczbaPunktow);

		/*double x[4] = { 0, 0.025, 0.025, 0 };
		double y[4] = { 0, 0, 0.025, 0.025 };*/
		double x[4];
		double y[4];

		for (int i = 0; i < 4; i++) {
			x[i] = nodes[ID[i] - 1].x;
			y[i] = nodes[ID[i] - 1].y;
		}

		int n = liczbaPunktow * liczbaPunktow;

		Pw punktyIWagi = metodaGaussa(2, 2);
		
		double** dNidx = new double* [n];
		double** dNidy = new double* [n];

		for (int i = 0; i < n; i++) {
			dNidx[i] = new double[4];
			dNidy[i] = new double[4];
		}

		for (int k = 0; k < n; k++) {
			double dxdEta = element.tabEta[k][0] * x[0] + element.tabEta[k][1] * x[1] + element.tabEta[k][2] * x[2] + element.tabEta[k][3] * x[3];
			double dxdKsi = element.tabKsi[k][0] * x[0] + element.tabKsi[k][1] * x[1] + element.tabKsi[k][2] * x[2] + element.tabKsi[k][3] * x[3];

			double dydEta = element.tabEta[k][0] * y[0] + element.tabEta[k][1] * y[1] + element.tabEta[k][2] * y[2] + element.tabEta[k][3] * y[3];
			double dydKsi = element.tabKsi[k][0] * y[0] + element.tabKsi[k][1] * y[1] + element.tabKsi[k][2] * y[2] + element.tabKsi[k][3] * y[3];

			double wyznacznik[2][2] = { {dxdKsi, dydKsi}, {dxdEta, dydEta} };


			double det = dxdKsi * dydEta - dydKsi * dxdEta;


			for (int j = 0; j < n; j++) {
				dNidx[k][j] = 1.0 / det * (wyznacznik[0][0] * element.tabKsi[k][j] + wyznacznik[0][1] * element.tabEta[k][j]);
				dNidy[k][j] = 1.0 / det * (wyznacznik[1][0] * element.tabKsi[k][j] + wyznacznik[1][1] * element.tabEta[k][j]);
			}


			/*for (int i = 0; i < n; i++) {
				for (int j = 0; j < 4; j++) {
					cout << dNidx[i][j] << " ";
				}
				cout << endl;
			}

			for (int i = 0; i < n; i++) {
				for (int j = 0; j < 4; j++) {
					cout << dNidy[i][j] << " ";
				}
				cout << endl;
			}*/

			for (int i = 0; i < n; i++) {
				for (int j = 0; j < 4; j++) {
					H[i][j] += conductivity * (dNidx[k][i] * dNidx[k][j] + dNidy[k][i] * dNidy[k][j]) * det;

				}
			}
		}

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				HBC[i][j] = 0.0;
				P[i] = 0.0;
			}
		}

		//BC
		for (int k = 0; k < 4; k++) {
			if (nodes[ID[k] - 1].BC == 1 && nodes[ID[(k + 1) % 4] - 1].BC == 1) {
				//cout << "HBC krawedz " << k + 1 << endl;

				double x1 = nodes[ID[k] - 1].x;
				double y1 = nodes[ID[k] - 1].y;
				double x2 = nodes[ID[(k + 1) % 4] - 1].x;
				double y2 = nodes[ID[(k + 1) % 4] - 1].y;

				double l = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
				double detJ = l * 0.5;
				
				//cout << "wyznacznik j" << detJ << endl;
				//tot = 1200;
				for (int z = 0; z < 4; z++) {
					for (int i = 0; i < 4; i++) {
						for (int j = 0; j < liczbaPunktowIntegracji; j++) {
							HBC[z][i] += alfa * punktyIWagi.w[j] * element.surface[k].N[j][z] * element.surface[k].N[j][i] * detJ;
						}
					}
				}

				for (int i = 0; i < liczbaPunktowIntegracji; i++) {
					for (int j = 0; j < 4; j++) {
						P[j] += alfa * (punktyIWagi.w[i] * element.surface[k].N[i][j] * tot) * detJ;
					}
				}
				
			}
		}

		for (int i = 0; i < n; i++)
		{
			delete[] dNidx[i];
			delete[] dNidy[i];

		}

		delete[] dNidx;
		delete[] dNidy;
	}

	~Element()
	{
		for (int i = 0; i < liczbaPunktowIntegracji * liczbaPunktowIntegracji; i++)
		{
			delete[] H[i];

		}

		delete[] H;

	}

	void showData(int index) {

		
		std::cout << "H" << index << endl;
		for (int i = 0; i < 4; i++)
			{
			for (int j = 0; j < 4; j++)
			{
				std::cout << H[i][j] << " ";
			}

			std::cout << "\n";
		}


		std::cout << endl << "P" << index << endl;
		for (int i = 0; i < 4; i++)
		{
			std::cout << P[i] << "\t";
		}

		std::cout << endl << "HBC" << index << endl;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				std::cout << HBC[i][j] << "\t";
			}

			std::cout << "\n";
		}
		std::cout << endl;
	}
};

#endif
