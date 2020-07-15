#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

using namespace std;

const int numWalkers = 100;  //NUMBER OF WALKERS
const int numSteps = 7;  //NUMBER OF STEPS

int main(int argc, char const* argv[])
{
    
    //SET OUTPUT TEXT FILES
    ofstream probOut;
    ofstream posMean;
    probOut.open("positionVSProbabilityOutput.csv");
    posMean.open("positionVSMeanSquareOutput.csv");

    cout << "NUMBER OF WALKERS: " << numWalkers << endl;
    cout << "NUMBER OF STEPS: " << numSteps << endl << endl;

    //aStepSquared IS WHERE A STEP OF A WALKER IS SQUARED 
    int* aStepSquared = new int[numSteps]; 
    for (int itime = 0; itime < numSteps; itime++) {
        aStepSquared[itime] = 0;
    }

    //SET PROBABILITY VARIABLES
    double* probability = new double[2 * numSteps + 1]; //ARRAY TO HOLD PROBABILITY OF EACH POSITION TO BE VISITED
    int* countNumberOfVisits = new int[numSteps * 2 + 1]; //ARRAY TO COUNT THE NUMBER OF VISITS PER POSITION

    for (int position = 0; position < (2 * numSteps + 1); position++) {
        countNumberOfVisits[position] = 0;
    }

    //FOR LOOP FOR EACH WALKER
    for (int kWalkr = 0; kWalkr < numWalkers; kWalkr++) {
        int currentPosition = 0; //CURRENT POSITION OF WALKER

        //FOR LOOP - FOR EACH STEP TAKEN
        for (int itime = 0; itime < numSteps; itime++)
        {
            //HEAD OR TAIL CONCEPT
            int randomNum = rand() % 2;

            //==============================
            if (randomNum == 1)
                currentPosition++;
            else 
                currentPosition--;
            //==============================
            cout << "==== New position is: " << currentPosition << endl;

            //COMPUTE aStepSquared[itime]
            aStepSquared[itime] = aStepSquared[itime] + (currentPosition * currentPosition);
            cout << "==== The squared displacement is x^2 = " << aStepSquared[itime] << endl << endl;
        }
        cout << endl;
        //INCREMENT countNumberOfVisits[] EVERYTIME currentPosition IS VISITED
        countNumberOfVisits[currentPosition + numSteps] += 1;
    }

    //Mean Square
    cout << endl;
    for (int itime = 0; itime < numSteps; itime++) {
        double meanSquare = (double)aStepSquared[itime] / (double)numWalkers;
        posMean << itime << ", " << meanSquare << endl;
    }

    //PROBABILITY
    cout << endl;
    cout << "Probability Values:" << endl << endl;
    for (int position = 0; position < (2 * numSteps + 1); position++)
    {
        //COMPUTE PROBABILITY OF POSITION i
        probability[position] = (double)countNumberOfVisits[position] / (double)numWalkers;
        if (countNumberOfVisits[position] != 0) {
            cout << "Total Steps and Probability of Walker To Be At Distance " << position - numSteps << ": " << countNumberOfVisits[position] << ", " << probability[position];
            probOut << position - numSteps << ", " << probability[position] << endl;  // Output(x,y): POSITION VS PROBABILITY
        }
    }
    probOut.close();
    posMean.close();
    return 0;
}
