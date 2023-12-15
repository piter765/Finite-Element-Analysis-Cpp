#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include "metodaGaussa.h"
#include "grid.h"
#include "elementUniwersalny.h"
#include "soe.h"

#define NUM_OF_POINTS 3
using namespace std;

int getDataFromFile(string filename, GlobalData& globalData, Grid& grid);
void generateParaViewFile(Grid grid, int i, Soe& soe);

int main() {
	std::cout << setprecision(10);
	GlobalData globalData = GlobalData();
	Grid grid = Grid();


	getDataFromFile("Test2_4_4_MixGrid.txt", globalData, grid);
	//showGlobalData(globalData);
	showGrid(grid);

	/*cout << metodaGaussa(1, 2) << endl;
	cout << metodaGaussa(2, 2) << endl;
	cout << metodaGaussa(2, 3) << endl;*/

	//ElementUniwersalny element = ElementUniwersalny(2);
	//element.printElementUniversal();

	//ElementUniwersalny elU = ElementUniwersalny(2);

	//elU.printElementUniversal();
	// 
	//Element el = Element(2);

	int iterations = globalData.simulationTime / globalData.simulationStepTime;

	Soe soe = Soe(grid, globalData);
	soe.calculate();

	generateParaViewFile(soe.grid, 0, soe); //dla initialTemp

	for (int i = 1; i <= iterations; i++) {
		soe.computeTempP();
		soe.calculateTemperature();
		/*if (i < 3) {
			soe.showSoe();
		}*/
		
		
		for (int k = 0; k < soe.grid.nE; k++) {
			std::cout << "temp!!!!: " << soe.tabNewTemperatures[k] << " ";
			
		}
	

		generateParaViewFile(soe.grid, i, soe);

	}


	return 0;
}

int getDataFromFile(string filename, GlobalData& globalData, Grid& grid) {
	ifstream inputFile(filename);
	if (!inputFile.is_open()) {
		cerr << "Failed to open the input file." << endl;
		return 0;
	}

	string line;
	int nodesNumber = 0;
	int elementsNumber = 0;

	while (getline(inputFile, line)) {
		istringstream iss(line);
		string keyword;
		iss >> keyword;
		if (keyword == "SimulationTime")
			iss >> globalData.simulationTime;
		else if (keyword == "SimulationStepTime")
			iss >> globalData.simulationStepTime;
		else if (keyword == "Conductivity")
			iss >> globalData.conductivity;
		else if (keyword == "Alfa")
			iss >> globalData.alfa;
		else if (keyword == "Tot")
			iss >> globalData.tot;
		else if (keyword == "InitialTemp")
			iss >> globalData.initialTemp;
		else if (keyword == "Density")
			iss >> globalData.density;
		else if (keyword == "SpecificHeat")
			iss >> globalData.specificHeat;
		else if (keyword == "Nodes") {
			string x;
			iss >> x >> grid.nN;
		}
		else if (keyword == "Elements") {
			string x;
			iss >> x >> grid.nE;
		}
		else if (keyword == "*Node") {
			string nvm;
			grid.nodes = new Node[grid.nN];
			for (int i = 0; i < grid.nN; i++) {
				string x, y;
				inputFile >> nvm >> x >> y;
				x.pop_back();
				y.pop_back();
				grid.nodes[i].x = stod(x);
				grid.nodes[i].y = stod(y);
				grid.nodes[i].t = globalData.initialTemp;
			}
		}
		else if (keyword == "*BC") {
			string x;

			while (inputFile >> x) {
				size_t commaPosition = x.find(",");
				if (commaPosition != std::string::npos) {
					x.erase(commaPosition, 1);
				}
				int index = stoi(x);
				grid.nodes[index - 1].BC = 1;
			}

			for (int i = 0; i < grid.nN; i++) {
				if (grid.nodes[i].BC != 1) {
					grid.nodes[i].BC = 0;
				}
				//cout << grid.nodes[i].BC << " ";
			}
		}

	}

	inputFile.close();
	ifstream inputFile2(filename);
	if (!inputFile2.is_open()) {
		cerr << "Failed to open the input file." << endl;
		return 0;
	}

	while (getline(inputFile2, line)) {
		istringstream iss(line);
		string keyword;
		iss >> keyword;
	    if (keyword == "*Element,") {
			string x;
			grid.elements = new Element[grid.nE];
			for (int i = 0; i < grid.nE; i++) {
				string one, two, three, four;
				inputFile2 >> x >> one >> two >> three >> four;
				one.pop_back();
				two.pop_back();
				three.pop_back();

				grid.elements[i].ID[0] = stod(one);
				grid.elements[i].ID[1] = stod(two);
				grid.elements[i].ID[2] = stod(three);
				grid.elements[i].ID[3] = stod(four);

				grid.elements[i].createElement(NUM_OF_POINTS, grid.nodes, globalData.conductivity, globalData.alfa, globalData.tot, globalData.density, globalData.specificHeat, globalData.simulationStepTime, globalData.initialTemp);
				//grid.elements[i].showData(i);
			}
		}
		
	}

}

void generateParaViewFile(Grid grid, int i, Soe& soe) {
	std::string fileName = "file_" + std::to_string(i) + ".vtk";

	ofstream outputFile(fileName);
	if (outputFile.is_open()) {

		/*# vtk DataFile Version 2.0
			Unstructured Grid Example
			ASCII
			DATASET UNSTRUCTURED_GRID*/
		outputFile << "# vtk DataFile Version 2.0\nUnstructured Grid Example\nASCII\nDATASET UNSTRUCTURED_GRID\n\n";
		outputFile << "POINTS " + std::to_string(grid.nN) + " float\n";

		for (int i = 0; i < grid.nN; i++) {
			outputFile << grid.nodes[i].x << " " << grid.nodes[i].y << " 0\n";
		}
		outputFile << "\n";

		outputFile << "CELLS " << grid.nE << " " << grid.nE * 5 << "\n";
		for (int i = 0; i < grid.nE; i++) {
			outputFile << "4 ";
			for (int j = 0; j < 4; j++) { // na razie 4 tyle co w ID jest
				outputFile << grid.elements[i].ID[j] - 1 << " ";
			}
			outputFile << "\n";
		}
		outputFile << "\n";

		outputFile << "CELL_TYPES " << grid.nE << "\n";
		for (int i = 0; i < grid.nE; i++) {
			outputFile << 9 << "\n";
		}
		outputFile << "\n";

		outputFile << "POINT_DATA " << grid.nN << "\n";
		outputFile << "SCALARS Temp float 1" << "\n";
		outputFile << "LOOKUP_TABLE default " << "\n";

		for (int i = 0; i < grid.nN; i++) {
			outputFile << soe.tabNewTemperatures[i] << "\n";
		}

		outputFile.close();
		std::cout << "Data has been written to the file." << std::endl;
	}
	else {
		std::cerr << "Unable to open file for writing." << std::endl;
	}
}