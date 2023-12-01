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

using namespace std;

int getDataFromFile(string filename, GlobalData& globalData, Grid& grid);

int main() {
	std::cout << setprecision(10);
	GlobalData globalData = GlobalData();
	Grid grid = Grid();


	getDataFromFile("Test1_4_4.txt", globalData, grid);
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

	Soe soe = Soe(grid, globalData);


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

				grid.elements[i].createElement(2, grid.nodes, globalData.conductivity, globalData.alfa, globalData.tot, globalData.density, globalData.specificHeat, globalData.simulationStepTime, globalData.initialTemp);
				grid.elements[i].showData(i);
			}
		}
		
	}

}

