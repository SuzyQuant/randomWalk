#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <string>
#include <vector>
using namespace std;


int main(int argc, char const* argv[])
{
    int numWalkers = 10000;
    int numSteps = 10;

    //Read binary values from binary file and store all values in an array except first column.
	//http://www.cplusplus.com/doc/tutorial/files/
    //http://courses.cs.vt.edu/cs2604/fall01/binio.html
    ifstream positions("positions.bin", ios::in | ios::binary);
    int* binPositions = new int[numWalkers * numSteps];
    int bp = 0;

    for (int row = 0; row < numSteps; row++) {
        for (int col = 0; col <= numWalkers; col++) {
            int val;
            positions.read((char*)&val, sizeof(int));   //extracted from positions then store it to val
                                                        //&val : Initial byte of an object stored in file. http://www.tutorialdost.com/Cpp-Programming-Tutorial/62-Cpp-File-Handling-IO-Read-Write-Object.aspx
                                                        //sizeof(int) : size of object represents the total number of bytes to be read from initial byte.
            //check if not first column
            if (col > 0) {
                binPositions[bp++] = val;
            }
        }
    }

    //Offstream files
    ofstream posMeanCSV;
    ofstream posMean;
    ofstream probOutCSV;
    ofstream probOut;

    posMeanCSV.open("quantitiesVStime.csv");
    posMean.open("quantitiesVStime.txt");
    probOutCSV.open("probabilityVSposition.csv");
    probOut.open("probabilityVSposition.txt");

    string stringHolder;

    //Variables for Squared Distance
    double q1 = 1 / sqrt(numSteps), q2 = 1, q3 = sqrt(numSteps);

    cout << "q1 " << q1 << endl;
    cout << "q2 " << q2 << endl;
    cout << "q3 " << q3 << endl;

    int* aStepSquared = new int[numSteps];
    int* aStepSquared2 = new int[numSteps];
    double* cosRq1 = new double[numSteps];
    double* cosRq2 = new double[numSteps];
    double* cosRq3 = new double[numSteps];
    double* cosSqrdRq1 = new double[numSteps];
    double* cosSqrdRq2 = new double[numSteps];
    double* cosSqrdRq3 = new double[numSteps];

    for (int itime = 0; itime < numSteps; itime++) {
        aStepSquared[itime] = 0;
        aStepSquared2[itime] = 0;
        cosRq1[itime] = 0.0;
        cosRq2[itime] = 0.0;
        cosRq3[itime] = 0.0;
        cosSqrdRq1[itime] = 0.0;
        cosSqrdRq2[itime] = 0.0;
        cosSqrdRq3[itime] = 0.0;
    }

    //SET PROBABILITY VARIABLES
    float* probability = new float[2 * numSteps + 1]; //ARRAY TO HOLD PROBABILITY OF EACH POSITION TO BE VISITED
    int* countNumberOfVisits = new int[numSteps * 2 + 1]; //ARRAY TO COUNT THE NUMBER OF VISITS PER POSITION

    for (int i = 0; i < (2 * numSteps + 1); i++) {
        countNumberOfVisits[i] = 0;
        probability[i] = 0.0;
    }

    //=====================================================//
    for (int itime = 0; itime < numSteps; itime++) {

        for (int iwalker = 0; iwalker < numWalkers; iwalker++) {

            int currentPosition = binPositions[itime * numWalkers + iwalker];

            //COMPUTE aStepSquared[itime] and cosRq1, cosRq2, cosRq3
            aStepSquared[itime] += (currentPosition * currentPosition);
            aStepSquared2[itime] += pow(currentPosition, 4);

            // cos (q*R)
            cosRq1[itime] += cos(q1 * currentPosition);
            cosRq2[itime] += cos(q2 * currentPosition);
            cosRq3[itime] += cos(q3 * currentPosition);

            //cos^2 (q*R)
            cosSqrdRq1[itime] += cos(q1 * currentPosition) * cos(q1 * currentPosition);
            cosSqrdRq2[itime] += cos(q2 * currentPosition) * cos(q2 * currentPosition);
            cosSqrdRq3[itime] += cos(q3 * currentPosition) * cos(q3 * currentPosition);
        }
        //Increase the currentPosition's number of visits
    }

    //PROBABILITY
    int itime = numSteps - 1; //final position only is necessary for probability
    for (int iwalker = 0; iwalker < numWalkers; iwalker++) {
            int currentPosition = binPositions[itime * numWalkers + iwalker];
            countNumberOfVisits[currentPosition + numSteps] += 1;
       
    }

    //Mean Square
    for (int itime = 0; itime < numSteps; itime++) {
        double meanSquare = (double)aStepSquared[itime] / (double)numWalkers;
        double meanSquare2 = (double)aStepSquared2[itime] / (double)numWalkers;
        //<cos(q*R)>
        cosRq1[itime] = (double)cosRq1[itime] / (double)numWalkers;
        cosRq2[itime] = (double)cosRq2[itime] / (double)numWalkers;
        cosRq3[itime] = (double)cosRq3[itime] / (double)numWalkers;
        //<cos^2(q*R)>
        cosSqrdRq1[itime] = (double)cosSqrdRq1[itime] / (double)numWalkers;
        cosSqrdRq2[itime] = (double)cosSqrdRq2[itime] / (double)numWalkers;
        cosSqrdRq3[itime] = (double)cosSqrdRq3[itime] / (double)numWalkers;
        //standard deviation of the mean
        double stdevMeanSquare = (double)sqrt(meanSquare2 - meanSquare * meanSquare);
        double stdErrorOfMean = (double)(stdevMeanSquare / sqrt(numWalkers));
        double stdevCosSqrdRq1 = (double)sqrt(cosSqrdRq1[itime] - (cosRq1[itime] * cosRq1[itime]));
        double stdErrorCosSqrdRq1 = (double)(stdevCosSqrdRq1 / sqrt(numWalkers));
        double stdevCosSqrdRq2 = (double)sqrt(cosSqrdRq2[itime] - (cosRq2[itime] * cosRq2[itime]));
        double stdErrorCosSqrdRq2 = (double)(stdevCosSqrdRq2 / sqrt(numWalkers));
        double stdevCosSqrdRq3 = (double)sqrt(cosSqrdRq3[itime] - (cosRq3[itime] * cosRq3[itime]));
        double stdErrorCosSqrdRq3 = (double)(stdevCosSqrdRq3 / sqrt(numWalkers));

        //long double exp1 = exp(-(itime + 1.0) / (2.0 * numSteps));
        //long double exp2 = exp(-(itime + 1.0) / 2.0);
        //long double exp3 = exp(-((itime + 1.0) * numSteps) / 2.0);
         
        long double exp1 = exp(-((itime + 1.0) / 2.0 ) * q1 * q1);
        long double exp2 = exp(-((itime + 1.0) / 2.0) * q2 * q2);
        long double exp3 = exp(-((itime + 1.0) / 2.0) * q3 * q3);

        cout << "exp" << endl;

        // .txt file
        posMean << itime << " " << meanSquare << " " << itime + 1 << " " << stdErrorOfMean << " " << cosRq1[itime] << " " << cosRq2[itime] << " " << cosRq3[itime] << " " << stdErrorCosSqrdRq1 << " " << stdErrorCosSqrdRq2 << " " << stdErrorCosSqrdRq3 << " " << setprecision(9) << exp1 << " " << exp2 << " " << exp3 << endl;
        //posMean << itime << " " << meanSquare << " " << itime + 1 << " " << stdErrorOfMean << " " << cosRq1[itime] << " " << cosRq2[itime] << " " << cosRq3[itime] << " " << stdErrorCosSqrdRq1 << " " << stdErrorCosSqrdRq2 << " " << stdErrorCosSqrdRq3 << endl;
        
        posMeanCSV << itime << "," << meanSquare << "," << itime + 1 << "," << stdErrorOfMean << "," << cosRq1[itime] << "," << cosRq2[itime] << "," << cosRq3[itime] << "," << stdErrorCosSqrdRq1 << "," << stdErrorCosSqrdRq2 << "," << stdErrorCosSqrdRq3 << ", " << setprecision(9) << exp1 << ", " << exp2 << ", " << exp3 << endl;

    }

    cout << "q1 " << q1 << endl;
    cout << "q2 " << q2 << endl;
    cout << "q3 " << q3 << endl;

    //OTHER PROBABILITY REQUESTED
    //UTILIZING ONLY THE VISITS ACTUALLY DONE
    //==================================
    for (int i = 0; i < (2 * numSteps + 1); i += 2)
    {
        probability[i] = (double)countNumberOfVisits[i] / (double)numWalkers;
        probOut << i - numSteps << " " << probability[i] << endl;  // Output(x,y): POSITION VS PROBABILITY
        probOutCSV << i - numSteps << ", " << probability[i] << endl;  // Output(x,y): POSITION VS PROBABILITY
    }

    posMeanCSV.close();
    posMean.close();
    probOut.close();
    probOutCSV.close();

    return 0;
}
