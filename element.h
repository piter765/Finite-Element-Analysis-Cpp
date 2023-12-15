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
	int liczbaPunktowIntegracji = 0;
	double P[4];
	double C[4][4];

	Element() : liczbaPunktowIntegracji() {
		for (int i = 0; i < 4; i++) {
			ID[i] = 0;
			for (int j = 0; j < 4; j++) {
				H[i][j] = 0.0;
				HBC[i][j] = 0.0;
				P[i] = 0.0;
				C[i][j] = 0;
			}
		}
	}

	void createElement(int liczbaPunktow, Node* nodes, double conductivity, double alfa, double tot, double density, double specificHeat, double simulationStepTime, double initialTemp) {

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

		Pw punktyIWagi = metodaGaussa(2, liczbaPunktowIntegracji);
		
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

			double wyznacznik[2][2] = { {dydEta, -dydKsi}, { -dxdEta, dxdKsi} };


			double det = dydEta * dxdKsi - (-dydKsi * -dxdEta);


			for (int j = 0; j < 4; j++) {
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

			/*double* tabWFor3point = {

			}*/
			double* tabW = new double[n]; //tablica wag 
			if (liczbaPunktow == 2) {
				tabW[0] = 1;//punktyIWagi.w[0] * punktyIWagi.w[0];
				tabW[1] = 1;// punktyIWagi.w[0] * punktyIWagi.w[1]; //potem
				tabW[2] = 1;// punktyIWagi.w[1] * punktyIWagi.w[1];

				tabW[3] = 1;// punktyIWagi.w[1] * punktyIWagi.w[2];
			}
			else if (liczbaPunktow == 3) {
				tabW[0] = punktyIWagi.w[0] * punktyIWagi.w[0];
				tabW[1] = punktyIWagi.w[0] * punktyIWagi.w[1];
				tabW[2] = punktyIWagi.w[0] * punktyIWagi.w[2];

				tabW[3] = punktyIWagi.w[1] * punktyIWagi.w[0];
				tabW[4] = punktyIWagi.w[1] * punktyIWagi.w[1];
				tabW[5] = punktyIWagi.w[1] * punktyIWagi.w[2];

				tabW[6] = punktyIWagi.w[2] * punktyIWagi.w[0];
				tabW[7] = punktyIWagi.w[2] * punktyIWagi.w[1];
				tabW[8] = punktyIWagi.w[2] * punktyIWagi.w[2];
			};

			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					H[i][j] += conductivity * tabW[k] * (dNidx[k][i] * dNidx[k][j] + dNidy[k][i] * dNidy[k][j]) * det;
					C[i][j] += density * specificHeat * tabW[k] * (element.tabN[k][i] * element.tabN[k][j]) * det;
				} // razy wagi pomnozyc jjeszcze 
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
		
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				H[i][j] += HBC[i][j];
				H[i][j] = H[i][j] + C[i][j] / simulationStepTime;
				
				//P[i] += C[i][j] / simulationStepTime * initialTemp;
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

		std::cout << "C" << index << endl;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				std::cout << C[i][j] << " ";
			}

			std::cout << "\n";
		}
		cout << endl;
	}
};

#endif
