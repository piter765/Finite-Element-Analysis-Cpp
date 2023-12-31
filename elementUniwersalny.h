#pragma once
#ifndef ELEMENT_UNIWERSALNY_H
#define ELEMENT_UNIWERSALNY_H

#include <cmath>
#include "metodaGaussa.h"

struct Surface {
    double** N;
};


struct ElementUniwersalny
{
    int punktyIntegracji;
    double** tabKsi = NULL;
    double** tabEta = NULL;
    double** tabN = NULL;

    Surface surface[4];

    ElementUniwersalny(int liczbaPunktow)
    {
        punktyIntegracji = liczbaPunktow;
        int punktyIntegracjiKwadrat = pow(punktyIntegracji, 2);

        Pw punktyIWagi = metodaGaussa(2, liczbaPunktow);


        tabKsi = new double* [punktyIntegracjiKwadrat];
        tabEta = new double* [punktyIntegracjiKwadrat];
        tabN = new double* [punktyIntegracjiKwadrat];

        for (int i = 0; i < punktyIntegracjiKwadrat; i++)
        {
            tabKsi[i] = new double[4];
            tabEta[i] = new double[4];
            tabN[i] = new double[4];
        }

        int* iE = new int[punktyIntegracjiKwadrat];
        int* iKsi = new int[punktyIntegracjiKwadrat];

        for (int i = 0; i < punktyIntegracjiKwadrat; i++) {


            if (punktyIntegracji < 2 && punktyIntegracji > 4) return;

            if (punktyIntegracji == 2) {
                iE[0] = 0; iE[1] = 0; iE[2] = 1; iE[3] = 1;
                iKsi[0] = 0; iKsi[1] = 1; iKsi[2] = 0; iKsi[3] = 1;
            }
            else if (punktyIntegracji == 3) {
                iE[0] = 0; iE[1] = 0; iE[2] = 0;
                iE[3] = 1; iE[4] = 1; iE[5] = 1;
                iE[6] = 2; iE[7] = 2; iE[8] = 2;

                iKsi[0] = 0; iKsi[1] = 1; iKsi[2] = 2;
                iKsi[3] = 0; iKsi[4] = 1; iKsi[5] = 2;
                iKsi[6] = 0; iKsi[7] = 1; iKsi[8] = 2;
            }
            else if (punktyIntegracji == 4) {
                iE[0] = 0; iE[1] = 0; iE[2] = 0; iE[3] = 0;
                iE[4] = 1; iE[5] = 1; iE[6] = 1; iE[7] = 1;
                iE[8] = 2; iE[9] = 2; iE[10] = 2; iE[11] = 2;
                iE[12] = 3; iE[13] = 3; iE[14] = 3; iE[15] = 3;

                iKsi[0] = 0; iKsi[1] = 1; iKsi[2] = 2; iKsi[3] = 3;
                iKsi[4] = 0; iKsi[5] = 1; iKsi[6] = 2; iKsi[7] = 3;
                iKsi[8] = 0; iKsi[9] = 1; iKsi[10] = 2; iKsi[11] = 3;
                iKsi[12] = 0; iKsi[13] = 1; iKsi[14] = 2; iKsi[15] = 3;
            }

            tabKsi[i][0] = -0.25 * (1.0 - punktyIWagi.p[iE[i]]);
            tabKsi[i][1] = 0.25 * (1.0 - punktyIWagi.p[iE[i]]);
            tabKsi[i][2] = 0.25 * (1.0 + punktyIWagi.p[iE[i]]);
            tabKsi[i][3] = -0.25 * (1.0 + punktyIWagi.p[iE[i]]);

            tabEta[i][0] = -0.25 * (1.0 - punktyIWagi.p[iKsi[i]]);
            tabEta[i][1] = -0.25 * (1.0 + punktyIWagi.p[iKsi[i]]);
            tabEta[i][2] = 0.25 * (1.0 + punktyIWagi.p[iKsi[i]]);
            tabEta[i][3] = 0.25 * (1.0 - punktyIWagi.p[iKsi[i]]);

        }
        //surface
        double bcNodes[4][2][2] = {
             {{-sqrt(1.0 / 3.0), -1.0}, {sqrt(1.0 / 3.0), -1.0}},
             {{1.0, -sqrt(1.0 / 3.0)}, {1.0, sqrt(1.0 / 3.0)}},
             {{sqrt(1.0 / 3.0), 1.0}, {-sqrt(1.0 / 3.0), 1.0}},
             {{-1.0, sqrt(1.0 / 3.0)}, {-1.0, -sqrt(1.0 / 3.0)}}
        };

        double bcNodes2[4][3][2] = {
            {{-sqrt(3.0 / 5.0), -1.0}, {0.0, -1.0}, {sqrt(3.0 / 5.0), -1.0}},
             {{1.0, -sqrt(3.0 / 5.0)}, {1.0, 0.0}, {1.0, sqrt(3.0 / 5.0)}},
             {{sqrt(3.0 / 5.0), 1.0}, {0.0, 1.0}, {-sqrt(3.0 / 5.0), 1.0}},
             {{-1.0, sqrt(3.0 / 5.0)}, {-1.0, 0.0}, {-1.0, -sqrt(3.0 / 5.0)}}
        };

        for (int i = 0; i < 4; i++) {
            surface[i].N = new double* [punktyIntegracji];

            for (int j = 0; j < punktyIntegracji; j++) {
                surface[i].N[j] = new double[4];
            
                if (punktyIntegracji == 2) {
                    surface[i].N[j][0] = (1.0 / 4.0) * (1 - bcNodes[i][j][0]) * (1 - bcNodes[i][j][1]);
                    surface[i].N[j][1] = (1.0 / 4.0) * (1 + bcNodes[i][j][0]) * (1 - bcNodes[i][j][1]);
                    surface[i].N[j][2] = (1.0 / 4.0) * (1 + bcNodes[i][j][0]) * (1 + bcNodes[i][j][1]);
                    surface[i].N[j][3] = (1.0 / 4.0) * (1 - bcNodes[i][j][0]) * (1 + bcNodes[i][j][1]);
                }
                else if (punktyIntegracji == 3) {
                    surface[i].N[j][0] = (1.0 / 4.0) * (1 - bcNodes2[i][j][0]) * (1 - bcNodes2[i][j][1]);
                    surface[i].N[j][1] = (1.0 / 4.0) * (1 + bcNodes2[i][j][0]) * (1 - bcNodes2[i][j][1]);
                    surface[i].N[j][2] = (1.0 / 4.0) * (1 + bcNodes2[i][j][0]) * (1 + bcNodes2[i][j][1]);
                    surface[i].N[j][3] = (1.0 / 4.0) * (1 - bcNodes2[i][j][0]) * (1 + bcNodes2[i][j][1]);
                }
               /* cout << surface[i].N[j][0] << " ";
                cout << surface[i].N[j][1] << " ";
                cout << surface[i].N[j][2] << " ";
                cout << surface[i].N[j][3] << " ";
                cout << endl;*/
            
            }
            cout << endl;
        }
        //

            delete[] iE;
            delete[] iKsi;


            double nodesC[4][2] = {
             {-1.0 / sqrt(3.0), -1.0 / sqrt(3.0)},
             {sqrt(1.0 / 3.0), -1.0 / sqrt(3.0)},
             {-1.0 / sqrt(3.0), sqrt(1.0 / 3.0)},
             {sqrt(1.0 / 3.0), sqrt(1.0 / 3.0)}
            };

            double nodesC2[9][2] = {
                {-sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0)}, {0.0, -sqrt(3.0 / 5.0)}, {sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0)},
                {-sqrt(3.0 / 5.0), 0.0}, {0.0, 0.0}, {sqrt(3.0 / 5.0), 0.0},
                {-sqrt(3.0 / 5.0), sqrt(3.0 / 5.0)}, {0.0, sqrt(3.0 / 5.0)}, {sqrt(3.0 / 5.0), sqrt(3.0 / 5.0)}
            };

            for (int i = 0; i < punktyIntegracjiKwadrat; i++) {
                if (punktyIntegracji == 2) {
                    tabN[i][0] = (1.0 / 4.0) * (1.0 - nodesC[i][0]) * (1.0 - nodesC[i][1]);
                    tabN[i][1] = (1.0 / 4.0) * (1.0 + nodesC[i][0]) * (1.0 - nodesC[i][1]);
                    tabN[i][2] = (1.0 / 4.0) * (1.0 + nodesC[i][0]) * (1.0 + nodesC[i][1]);
                    tabN[i][3] = (1.0 / 4.0) * (1.0 - nodesC[i][0]) * (1.0 + nodesC[i][1]);
                }
                else if (punktyIntegracji == 3) {
                    tabN[i][0] = (1.0 / 4.0) * (1.0 - nodesC2[i][0]) * (1.0 - nodesC2[i][1]);
                    tabN[i][1] = (1.0 / 4.0) * (1.0 + nodesC2[i][0]) * (1.0 - nodesC2[i][1]);
                    tabN[i][2] = (1.0 / 4.0) * (1.0 + nodesC2[i][0]) * (1.0 + nodesC2[i][1]);
                    tabN[i][3] = (1.0 / 4.0) * (1.0 - nodesC2[i][0]) * (1.0 + nodesC2[i][1]);
                }
            }
    }


    void printElementUniversal()
    {
        std::cout << "pochodna Ksi:\n";
        for (int i = 0; i < pow(punktyIntegracji, 2); i++)
        {
            for (int j = 0; j < 4; j++)
            {
                std::cout << tabKsi[i][j] << " ";
            }

            std::cout << "\n";
        }

        std::cout << "pochodna Eta:\n";
        for (int i = 0; i < pow(punktyIntegracji,2); i++)
        {
            for (int j = 0; j < 4; j++)
            {
                std::cout << tabEta[i][j] << " ";
            }

            std::cout << "\n";
        }
    }

    ~ElementUniwersalny()
    {
        for (int i = 0; i < pow(punktyIntegracji, 2); i++)
        {
            delete[] tabKsi[i];
            delete[] tabEta[i];
        }

        delete[] tabKsi;
        delete[] tabEta;
    }

};



#endif