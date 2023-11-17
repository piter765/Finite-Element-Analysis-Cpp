#pragma once

#ifndef GRID_H
#define GRID_H

#include "elementUniwersalny.h"
#include "element.h"


struct GlobalData {
	double simulationTime;
	double simulationStepTime;
	double conductivity;
	double alfa;
	double tot;
	double initialTemp;
	double density;
	double specificHeat;
};

struct Grid {
	int nN = 0;
	int nE = 0;
	Node* nodes;
	Element* elements;
};


void showGlobalData(GlobalData globalData) {
	cout << "SimulationTime: " << globalData.simulationTime << endl;
	cout << "SimulationStepTime: " << globalData.simulationStepTime << endl;
	cout << "Conductivity: " << globalData.conductivity << endl;
	cout << "Alfa: " << globalData.alfa << endl;
	cout << "Tot: " << globalData.tot << endl;
	cout << "InitialTemp: " << globalData.initialTemp << endl;
	cout << "Density: " << globalData.density << endl;
	cout << "SpecificHeat: " << globalData.specificHeat << endl;
}

void showGrid(Grid grid) {
	cout << "Nodes: " << endl;
	for (int i = 0; i < grid.nN; i++) {
		cout << "x: " << grid.nodes[i].x << " " << "y: " << grid.nodes[i].y << endl;
	}
	cout << endl;
	cout << "Elements: " << endl;
	for (int i = 0; i < grid.nE; i++) {
		cout << grid.elements[i].ID[0] << " " << grid.elements[i].ID[1] << " " << grid.elements[i].ID[2] << " " << grid.elements[i].ID[3] << endl;
	}
}

#endif

