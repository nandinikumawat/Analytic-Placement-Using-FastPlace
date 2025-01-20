// written by Suraj Dalvi, ECE Dept., Univ. of Minnesota
// Heavily modified by Kia Bazargan (renaming variables, adding
// comments, creating a clean version of the hyperedge data structure)


# include<iostream>
# include<stdio.h>
# include<map>
# include<vector>
#include <string.h>
#include <algorithm>
#include <cmath>
#include <fstream>


using namespace std;

struct ltstr
{
	bool operator()(const char* s1, const char* s2) const
	{
		return strcmp(s1, s2) < 0;
	}
};

map <const char*, int, ltstr> nodeNameToNodeNum_map;

vector<vector<double> > matQ;
int matrixSize = 0 ;
vector<double> vectorDX;
vector<double> vectorDY;

int numCellPins, 		// number of all terminals connected to the end points of (hyper) edges
numhyper, 			// number of edges and hyperedges
numCellsAndPads, 	// total number of movable cells (generall with names starting with a), and I/O pads (generally starting with p)
numCells_noPads;	// total number of movable cells

int* cellPinArray;		// array of cell names used as endpoints of (hyper) edges
// The size of the array is numCellPins
//vector<int> movableAndStarCells; // Vector to store only the movable and star cells
int* hEdge_idxToFirstEntryInPinArray;
int* hyperwts;			// (hyper) edge weights. The size of the array is numhyper
int* vertexSize;		// cell and I/O pad sizes. The size of the array is numCellsAndPads
vector<double> weightsNew; // contains all the changed weights after star and clique models

// vectorDX.resize(numCellsAndPads, 0.0);
// vectorDY.resize(numCellsAndPads, 0.0);
// weightsNew.resize(numhyper, 0.0);

struct SPinLocation {
	int x, y;
};


SPinLocation* pinLocations;



// Helper function for matrix-vector multiplication
vector<double> multiplyingMatrices(const vector<vector<double> >& mat, const vector<double>& vec) {
    vector<double> result(mat.size(), 0.0);
    for (int i = 0; i < mat.size(); ++i) {
        for (int j = 0; j < mat[i].size(); ++j) {
            result[i] += mat[i][j] * vec[j];
        }
    }
    return result;
}

double multiplicationValues(const vector<double>& p, const vector<double>& q) {
    double resultant = 0.0;
    for (int i = 0; i < p.size(); ++i) {
        resultant += p[i] * q[i];
    }
    return resultant;
}


// Conjugate Gradient Solver
vector<double> conjugateGradient(const vector<vector<double> >& Q, const vector<double>& b, double tol = 1e-6, int iterMax = 1000) {
    int n = b.size();
    vector<double> x(n, 0.0); // Initial guess
    vector<double> r = b;    // Initial residual (b - Qx, where x = 0)
    vector<double> p = r;    // Initial search direction
    double oldResultant = multiplicationValues(r, r);

    for (int i = 0; i < iterMax; ++i) {
        vector<double> K = multiplyingMatrices(Q, p);
        double alphaValue = oldResultant / multiplicationValues(p, K);
        for (int j = 0; j < n; ++j) {
            x[j] += alphaValue * p[j];
            r[j] -= alphaValue * K[j];
        }

    //         int n = b.size();
    // vector<double> x(n, 0.0); // Initial guess
    // vector<double> r = b;    // Initial residual (b - Qx, where x = 0)
    // vector<double> p = r;    // Initial search direction
    // double oldResultant = multiplicationValues(r, r);

		//         double newResultant = multiplicationValues(r, r);
        // if (sqrt(newResultant) < tol) {
        //     cout << "Converged in " << i + 1 << " iterations." << endl;
        //     break;
        // }
        // double valBeta = newResultant / oldResultant;
        // for (int j = 0; j < n; ++j) {
        //     p[j] = r[j] + valBeta * p[j];
        // }
        // oldResultant = newResultant;

        double newResultant = multiplicationValues(r, r);
        if (sqrt(newResultant) < tol) {
            cout << "Converged in " << i + 1 << " iterations." << endl;
            break;
        }
        double valBeta = newResultant / oldResultant;
        for (int j = 0; j < n; ++j) {
            p[j] = r[j] + valBeta * p[j];
        }
        oldResultant = newResultant;
    }
	return x;
}

double weightedWirelengthCalculation(const vector<vector<double> > &Q, const vector<double> &x, const vector<double> &y) {
    int n = Q.size();
    double wireLengthTotal = 0.0;

    // Iterate through each combination of nodes
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) { // Check only pairs (i, j) where i < j
            if (Q[i][j] < 0) { // Edge exists if Q[i][j] < 0
                double weight = -Q[i][j]; // Take the absolute weight

	// 			// Iterate through each combination of nodes
    // for (int i = 0; i < n; ++i) {
    //     for (int j = i + 1; j < n; ++j) { // Check only pairs (i, j) where i < j
    //         if (Q[i][j] < 0) { // Edge exists if Q[i][j] < 0
    //             double weight = -Q[i][j]; // Take the absolute weight
    //             double squaredDistanceCalculatiion = pow((-x[i]) - (-x[j]), 2) + pow((-y[i]) - (-y[j]), 2);
    //             wireLengthTotal += weight * squaredDistanceCalculatiion;
    //         }
    //     }
    // }
                double squaredDistanceCalculatiion = pow((-x[i]) - (-x[j]), 2) + pow((-y[i]) - (-y[j]), 2);
                wireLengthTotal += weight * squaredDistanceCalculatiion;
            }
        }
    }
    // Take the square root of the total sum
    return sqrt(wireLengthTotal);
}


// Function to perform cell spreading
void cellSpreading(vector<double>& xCoordinates, vector<double>& yCoordinates, 
                   double minValX, double maxValX, double minValY, double maxValY, int gridSize) {
    int numberOfXBins = gridSize;
    int numberOfYBins = gridSize;
    double widthOfBin = (maxValX - minValX) / numberOfXBins;
    double heightOfBin = (maxValY - minValY) / numberOfYBins;

    vector<vector<double>> binUtilization(numberOfXBins, vector<double>(numberOfYBins, 0.0));

    // Calculate bin utilization
    for (int i = 0; i < xCoordinates.size(); ++i) {
        int binX = min(max(static_cast<int>((xCoordinates[i] - minValX) / widthOfBin), 0), numberOfXBins - 1);
        int binY = min(max(static_cast<int>((yCoordinates[i] - minValY) / heightOfBin), 0), numberOfYBins - 1);
        binUtilization[binX][binY] += 1.0; // Assuming unit cell size for simplicity
    }

    // Cell shifting
    for (int i = 0; i < xCoordinates.size(); ++i) {
        int binX = min(max(static_cast<int>((xCoordinates[i] - minValX) / widthOfBin), 0), numberOfXBins - 1);
        int binY = min(max(static_cast<int>((yCoordinates[i] - minValY) / heightOfBin), 0), numberOfYBins - 1);

        double shiftX = 0.0, shiftY = 0.0;
        
        // Calculate shift based on neighboring bin utilization
        if (binX > 0) shiftX -= binUtilization[binX-1][binY];
        if (binX < numberOfXBins-1) shiftX += binUtilization[binX+1][binY];
        if (binY > 0) shiftY -= binUtilization[binX][binY-1];
        if (binY < numberOfYBins-1) shiftY += binUtilization[binX][binY+1];

        // Apply movement control
        double alpha = 0.5; // Adjust this value as needed
        xCoordinates[i] += alpha * shiftX * widthOfBin;
        yCoordinates[i] += alpha * shiftY * heightOfBin;
    }

    // Add spreading forces (simplified version)
    for (int i = 0; i < xCoordinates.size(); ++i) {
        double forceX = (xCoordinates[i] - (minValX + maxValX) / 2) * 0.1;
        double forceY = (yCoordinates[i] - (minValY + maxValY) / 2) * 0.1;
        xCoordinates[i] -= forceX;
        yCoordinates[i] -= forceY;
    }
}

double calculateWirelengthAfterSpreading(const vector<vector<double>>& matQ,
                                         const vector<double>& xCoordinates,
                                         const vector<double>& yCoordinates) {
    double wirelength = 0.0;

    // Iterate through the Q-matrix
    for (int i = 0; i < matQ.size(); ++i) {
        for (int j = i + 1; j < matQ[i].size(); ++j) {
            if (matQ[i][j] < 0) { // Check if there's an edge
                double dx = xCoordinates[i] - xCoordinates[j];
                double dy = yCoordinates[i] - yCoordinates[j];
                wirelength += -matQ[i][j] * sqrt(dx * dx + dy * dy); // Using negative Q for weights

            }
        }
    }
    return wirelength;
}

void writeCoordinatesToFile(const string& filename, 
                            const vector<double>& xCoordinates, 
                            const vector<double>& yCoordinates) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << " to write coordinates.\n";
        return;
    }

    for (int i = 0; i < xCoordinates.size(); ++i) {
        file << "Cell_" << i << " " << xCoordinates[i] << " " << yCoordinates[i] << "\n";
    }

    file.close();
    cout << "Coordinates written to " << filename << endl;
}

int parseIbmFile(char* inareFileName, char* innetFileName, char* inPadLocationFileName)
{
	char line[80];
	char vertexname[10], nodedesc[10];
	char hyperNet;
	int hyperweight;
	FILE* innetFile, * inareFile, * inPadLocationFile;

	innetFile = fopen(innetFileName, "r");
	// if (!innetFile) {
	// 	printf("ERROR: Cannot open input file %s.\n", innetFileName);
	// 	return -1;
	// }
	inareFile = fopen(inareFileName, "r");
	// if (!inareFile) {
	// 	printf("ERROR: Cannot open input file %s.\n", inareFileName);
	// 	return -1;
	// }
	inPadLocationFile = fopen(inPadLocationFileName, "r");
	// if (!inPadLocationFile) {
	// 	printf("ERROR: Cannot open input file %s.\n", inPadLocationFileName);
	// 	return -1;
	// }

	fscanf(innetFile, "%*d %d %d %d %d\n", &numCellPins, &numhyper, &numCellsAndPads, &numCells_noPads);
	
	numCells_noPads = numCells_noPads + 1;
	matQ.clear();
 	matQ.resize(numCellsAndPads, vector<double>(numCellsAndPads, 0.0));
	vectorDX.assign(numCellsAndPads, 0.0);
	vectorDY.assign(numCellsAndPads, 0.0);

	// 	fscanf(innetFile, "%*d %d %d %d %d\n", &numCellPins, &numhyper, &numCellsAndPads, &numCells_noPads);
	
	// numCells_noPads = numCells_noPads + 1;
	// matQ.clear();
 	// matQ.resize(numCellsAndPads, vector<double>(numCellsAndPads, 0.0));
	// vectorDX.assign(numCellsAndPads, 0.0);
	// vectorDY.assign(numCellsAndPads, 0.0);


	cout << "numCellPins, numhyper, numCellsAndPads, numCells_noPads = " << numCellPins << ", " << numhyper << ", " << numCellsAndPads << ", " << numCells_noPads << endl;

	vertexSize = (int*)malloc(numCellsAndPads * sizeof(int));

	for (int i = 0; i < numCellsAndPads; i++)
	{
		fscanf(inareFile, "%s %d", vertexname, &vertexSize[i]);
		//printf("\n%s %d",vertexname,vertexSize[i]);
		nodeNameToNodeNum_map[strdup(vertexname)] = i;
		// if cell name starts with "a", it's a movable cell
		// and if starts with "p", it is an I/O pad.
	}
	// Read (hyper)edges
	hEdge_idxToFirstEntryInPinArray = (int*)malloc((numhyper + 1) * sizeof(int)); // allocating memory for the array
	cellPinArray = (int*)malloc(numCellPins * sizeof(int));
	hyperwts = (int*)malloc(numhyper * sizeof(int));
	vectorDX.reserve(numCellsAndPads);
    vectorDY.reserve(numCellsAndPads);
    weightsNew.reserve(numhyper);

	if (hEdge_idxToFirstEntryInPinArray == NULL || cellPinArray == NULL || hyperwts == NULL)
	{
		printf("\nUnable to allocate memory");
		return -1;
	}
	int hypercount = 0;
	int pinCount = 0;

	while (fgets(line, 80, innetFile) != NULL)
	{
		sscanf(line, "%s %c %d", nodedesc, &hyperNet, &hyperweight);
		if (hyperNet == 's')
		{
			hEdge_idxToFirstEntryInPinArray[hypercount] = pinCount;
			hyperwts[hypercount] = hyperweight;
			++hypercount;
		}
		else
		{
			sscanf(line, "%s %d", nodedesc, &hyperweight);
		}
		cellPinArray[pinCount] = nodeNameToNodeNum_map[nodedesc];
		++pinCount;
	}
	hEdge_idxToFirstEntryInPinArray[hypercount] = numCellPins;

	// Read I/O Pad locations
	int numPads = numCellsAndPads - numCells_noPads;
	pinLocations = new SPinLocation[numPads];
	for (int i = 0; i < numPads; i++)
	{
		char padName[80];
		int x, y;
		fscanf(inPadLocationFile, "%s %d %d", padName, &x, &y);
		int padCellIdx = nodeNameToNodeNum_map[padName];
		pinLocations[padCellIdx - numCells_noPads].x = x;		// pay attention how it's indexed 
		pinLocations[padCellIdx - numCells_noPads].y = y;
		// Add x and y values to the vectors
		//vectorDX.push_back(static_cast<double>(x));
		//vectorDY.push_back(static_cast<double>(y));
	}


	matQ.resize(numCellsAndPads, vector<double>(numCellsAndPads, 0.0));

	for (int i = 0; i < numhyper; ++i)
	{
		weightsNew.push_back(static_cast<double>(hyperwts[i]));
	}


	map<int, int> countEle;

	// Iterate through cellPinArray
	for (int i = 0; i < numCellPins; ++i) {
		int currEle = cellPinArray[i];

		// Check if the element is already in the map
		auto it = countEle.find(currEle);
		if (it != countEle.end()) {
			// Increment the count if the element is already in the map
			it->second++;
		}
		else {
			// Insert the element with count 1 if it's not in the map
			countEle[currEle] = 1;

	// 		// Check if the element is already in the map
	// 	auto it = countEle.find(currEle);
	// 	if (it != countEle.end()) {
	// 		// Increment the count if the element is already in the map
	// 		it->second++;
	// 	}
	// 	else {
	// 		// Insert the element with count 1 if it's not in the map
	// 		countEle[currEle] = 1;
	// 	}
	// }
		}
	}
	for (const auto& entry : countEle) {
		//if (entry.second > 1) {
		//cout << "Element " << entry.first << " is repeated " << entry.second << " times." << endl;
		matQ[entry.first][entry.first] = static_cast<double>(entry.second);
		//}
	}

	//Initial Q matrix 
	matrixSize = numCells_noPads;
	


		for (int i = 0; i < numhyper; ++i) {
			// Get the starting and ending indices for the current hyperedge
			int startIndex = hEdge_idxToFirstEntryInPinArray[i];
			int endIndex = hEdge_idxToFirstEntryInPinArray[i + 1];

			// Check if startIndex and endIndex differ by 2
			if (endIndex - startIndex == 2) {

				// Access values in cellPinArray between startIndex and endIndex
				///cout << "Values in cellPinArray for hyperedge " << i << ": ";
				for (int j = startIndex; j < endIndex; j = j+2) {
					int indexC1 = cellPinArray[j];
					int indexC2 = cellPinArray[j+1];
					//cout << "indexC1 " << indexC1 << " ";
					//cout << "dog";

					// Check if cellIndex is less than numCells_noPads
					if (indexC1 < matrixSize && indexC2 < matrixSize) {
						// Update the matQ[cellIndex][cellIndex + 1] (assuming +1 is valid)
						matQ[indexC1][indexC2] = -1;
						matQ[indexC2][indexC1] = -1;

						break;
					}

				}
			}
			// cout << endl;
		}	
	//int y = numCells_noPads;

	int c = 1;
	int k;
	double weightValNew ;
	double gamma;
	double weightValOld;
	for (int i = 0; i < numhyper; i++)
	{
		if (hEdge_idxToFirstEntryInPinArray[i + 1] - hEdge_idxToFirstEntryInPinArray[i] == 3 )
		{
			//clique 
			int a = hEdge_idxToFirstEntryInPinArray[i + 1];
			int b = hEdge_idxToFirstEntryInPinArray[i];
			k = hEdge_idxToFirstEntryInPinArray[i + 1] - hEdge_idxToFirstEntryInPinArray[i];
			gamma = 1.0 / (k - 1);
			weightValOld = weightsNew[cellPinArray[hEdge_idxToFirstEntryInPinArray[i]]];
			weightValNew = gamma * weightValOld;
			weightsNew[cellPinArray[hEdge_idxToFirstEntryInPinArray[i]]] = weightValNew;
			//cout << "b = " << b << "a = " << a;
			//cout << cellPinArray[b] << " " << cellPinArray[b+1] << " " << cellPinArray[b+2];

			// //clique 
			// int a = hEdge_idxToFirstEntryInPinArray[i + 1];
			// int b = hEdge_idxToFirstEntryInPinArray[i];
			// k = hEdge_idxToFirstEntryInPinArray[i + 1] - hEdge_idxToFirstEntryInPinArray[i];
			// gamma = 1.0 / (k - 1);
			// weightValOld = weightsNew[cellPinArray[hEdge_idxToFirstEntryInPinArray[i]]];
			// weightValNew = gamma * weightValOld;
			// weightsNew[cellPinArray[hEdge_idxToFirstEntryInPinArray[i]]] = weightValNew;

					matQ[cellPinArray[b+1]][cellPinArray[b]] = matQ[cellPinArray[b+1]][cellPinArray[b]] - weightValNew;
					matQ[cellPinArray[b+2]][cellPinArray[b]] = matQ[cellPinArray[b+2]][cellPinArray[b]] - weightValNew;
					matQ[cellPinArray[b]][cellPinArray[b+1]] = matQ[cellPinArray[b]][cellPinArray[b+1]] -  weightValNew;

					matQ[cellPinArray[b+2]][cellPinArray[b + 1]] = matQ[cellPinArray[b+2]][cellPinArray[b + 1]] - weightValNew;
					matQ[cellPinArray[b]][cellPinArray[b+2]] = matQ[cellPinArray[b]][cellPinArray[b+2]] - weightValNew;
					matQ[cellPinArray[b+1]][cellPinArray[b+2]] = matQ[cellPinArray[b+1]][cellPinArray[b+2]] - weightValNew;
		}
		else if (hEdge_idxToFirstEntryInPinArray[i + 1] - hEdge_idxToFirstEntryInPinArray[i] > 3)
		{

			int a = hEdge_idxToFirstEntryInPinArray[i + 1];
			int b = hEdge_idxToFirstEntryInPinArray[i];
			//star model
			
			matrixSize = matrixSize + 1;
			//movableAndStarCells.push_back("as" + to_string(y));
			

			k = hEdge_idxToFirstEntryInPinArray[i + 1] - hEdge_idxToFirstEntryInPinArray[i];
			gamma = 1.0 / (k - 1);
			weightValOld = weightsNew[cellPinArray[hEdge_idxToFirstEntryInPinArray[i]]];
			weightValNew = k * gamma * weightValOld; // value that should go into the matQ
			weightsNew[cellPinArray[hEdge_idxToFirstEntryInPinArray[i]]] = weightValNew;
			for (int j = b; j < a; ++j)
			{
				if (cellPinArray[j] < numCells_noPads)
				{
					matQ[cellPinArray[j]][matrixSize - 1] = matQ[cellPinArray[j]][matrixSize - 1] - weightValNew;
					matQ[matrixSize - 1][cellPinArray[j]]= matQ[matrixSize - 1][cellPinArray[j]] - weightValNew;
					matQ[cellPinArray[j]][cellPinArray[j]] = matQ[cellPinArray[j]][cellPinArray[j]] - 1 + weightValNew; 
				}				
			}
			matQ[matrixSize - 1][matrixSize - 1] = weightValNew * (a-b);
			//cout << endl;
		}
	}
	cout << "matQ:" << endl;
	for (int i = 0; i < matrixSize; i++) {
		for (int j = 0; j < matrixSize; j++) {
			cout << matQ[i][j] << "\t";
		}
		cout << endl;
	}
	fclose(inPadLocationFile);
	fclose(innetFile);
	fclose(inareFile);

/****** Calculating dx and dy values ********/

	int X = 0;
	int Y = 0;
	double deg = 0;
	int countStarNode = 0;
	int nodeIdxVal = 0;
	int hypecount = 0;

	// vectorDX and vectorDY initialization 
	vector<double> vectorDX(matrixSize, 0.0);
	vector<double> vectorDY(matrixSize, 0.0);

	// Iterate through hyperedges
	for (int hyperedgeIndex = 0; hyperedgeIndex < numhyper; ++hyperedgeIndex) {
		int startPinIndex = hEdge_idxToFirstEntryInPinArray[hyperedgeIndex];
		int endPinIndex = hEdge_idxToFirstEntryInPinArray[hyperedgeIndex + 1];
		deg = endPinIndex - startPinIndex;

		// Vector to store node indices
		vector<int> nodeIndices;
		// Iterate through pins in the hyperedge
		for (int pinIndex = startPinIndex; pinIndex < endPinIndex; ++pinIndex) {
			nodeIdxVal = cellPinArray[pinIndex];
			nodeIndices.push_back(nodeIdxVal);
		}
		// Sort node indices in descending order
		sort(nodeIndices.rbegin(), nodeIndices.rend());

		// Update dx and dy based on the sorted order
		for (int i : nodeIndices) {
			// Check if the node corresponds to an I/O pads
			if (i >= numCells_noPads) {
				// Adjust index to match the pinLocations array
				int padIndex = i - numCells_noPads;
				X += pinLocations[padIndex].x;
				Y += pinLocations[padIndex].y;
			}
			else {
				// Update only the star node X and Y values
				if (deg > 3) {
					int indexStarNodeVal = numCells_noPads + countStarNode;

					if (indexStarNodeVal < matrixSize) {
						vectorDX[indexStarNodeVal] = -X * (deg / (deg - 1));
						vectorDY[indexStarNodeVal] = -Y * (deg / (deg - 1));
						hypecount++;
					}
				}
				else {
					// Update dx and dy based on the pad's location
					if (i < matrixSize) {
						vectorDX[i] -= X;
						vectorDY[i] -= Y;
					}
				}
			}

			// Increment countStarNode after processing all nodes in the hyperedge
			if (hypecount == nodeIndices.size() - 1) {
				countStarNode++;

				hypecount = 0;
			}
		}

		X = 0;
		Y = 0;
	}
	// Print the contents of dx
	cout << "dx Vector:" << endl;
	for (int i = 0; i < matrixSize; ++i) {
		cout << "Node " << i << ": " << vectorDX[i] << endl;
	}

	// Print the contents of dy
	cout << "dy Vector:" << endl;
	for (int i = 0; i < matrixSize; ++i) {
		cout << "Node " << i << ": " << vectorDY[i] << endl;
	}
	//cout << "matrix size " << matrixSize << endl;

	/***** Conjugate Gradient ********/


// Helper function for vector dot product
vector<double> x_coord = conjugateGradient(matQ,vectorDX);
for(int i=0; i<x_coord.size(); i++){
	cout << "x is: " << x_coord[i];
	cout << endl;
}
vector<double> y_coord = conjugateGradient(matQ,vectorDY);
for(int i=0; i<y_coord.size(); i++){
	cout << "y is: " << y_coord[i];
	cout << endl;
}

vector<double> xCoordinates = conjugateGradient(matQ, vectorDX);
vector<double> yCoordinates = conjugateGradient(matQ, vectorDY);

// Define chip area bounds (example values)
double minValX = *min_element(xCoordinates.begin(), xCoordinates.end());
double maxValX = *max_element(xCoordinates.begin(), xCoordinates.end());
double minValY = *min_element(yCoordinates.begin(), yCoordinates.end());
double maxValY = *max_element(yCoordinates.begin(), yCoordinates.end());

cout << "Total wirelength is: " << weightedWirelengthCalculation(matQ,x_coord,y_coord);
cout << endl;
writeCoordinatesToFile("coordinates_before_spreading.kiaPad", xCoordinates, yCoordinates);
cout << "Cell coordinates before spreading:\n";
for (int i = 0; i < xCoordinates.size(); ++i) {
    cout << "Cell " << i << ": (" << xCoordinates[i] << ", " << yCoordinates[i] << ")\n";
}
// Perform cell spreading
cellSpreading(xCoordinates, yCoordinates, minValX, maxValX, minValY, maxValY, 10); // Use a 10x10 grid

writeCoordinatesToFile("coordinates_after_spreading.kiaPad", xCoordinates, yCoordinates);
cout << "Cell coordinates after spreading:\n";
for (int i = 0; i < xCoordinates.size(); ++i) {
    cout << "Cell " << i << ": (" << xCoordinates[i] << ", " << yCoordinates[i] << ")\n";
}

double wirelengthAfterSpreading = calculateWirelengthAfterSpreading(matQ, xCoordinates, yCoordinates);
cout << "Wirelength after spreading: " << wirelengthAfterSpreading << endl;


return 0;
}


