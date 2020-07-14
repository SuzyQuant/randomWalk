#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

using namespace std;

const int numWalkers = 10;  //NUMBER OF WALKERS
const int numSteps = 50;  //NUMBER OF STEPS

int main(int argc, char const* argv[])
{
    string line;

    //SET OUTPUT TEXT FILES
    ofstream probOut;
    ofstream posNSteps;
    ofstream posMean;
    probOut.open("positionVSProbabilityOutput.csv");
    posNSteps.open("positionVSNumStepsOutput.csv");
    posMean.open("positionVSMeanSquareOutput.csv");

    cout << "NUMBER OF WALKERS: " << numWalkers << endl;
    cout << "NUMBER OF STEPS: " << numSteps << endl << endl;

    //ARRAY TO HOLD ALL X2AVE VALUES
    int* x2ave = new int[numSteps]; 
    for (int i = 0; i < numSteps; i++) {
        x2ave[i] = 0;
    }

    //SET PROBABILITY VARIABLES
    double* prb = new double[2 * numSteps + 1]; //ARRAY TO HOLD PROBABILITY OF EACH POSITION TO BE VISITED
    int* countNumberOfVisits = new int[numSteps * 2 + 1]; //ARRAY TO COUNT THE NUMBER OF VISITS PER POSITION
    for (int i = 0; i < (2 * numSteps + 1); i++) {
        countNumberOfVisits[i] = 0;
        prb[i] = 0.0;
    }

    //FOR LOOP FOR EACH WALKER
    for (int k = 0; k < numWalkers; k++) {
        int currentPosition = 0; //CURRENT POSITION OF WALKER

        //FOR LOOP - FOR EACH STEP TAKEN
        for (int i = 0; i < numSteps; i++)
        {
            //TOSS COIN IF HEAD OR TAIL
            int randomNum = rand() % 2;

            if (randomNum == 1) //IF HEADS
            {
                cout << "Walker " << k + 1 << " Step " << i + 1 << ":" << endl;
                currentPosition += 1;
                cout << "==== TOSS COIN SHOWS HEAD: x added so new position x is " << currentPosition << endl;
            }
            else //IF TAILS
            {
                cout << "Walker " << k + 1 << " Step " << i + 1 << ":" << endl;
                currentPosition += -1;
                cout << "==== COIN SHOWS TAIL: x subtracted so new position x is " << currentPosition << endl;
            }

            //INCREMENT countNumberOfVisits[i] EVERYTIME POSITION i IS VISITED
            countNumberOfVisits[currentPosition + numSteps] += 1;
            cout << "==== New position is: " << currentPosition << endl;

            //COMPUTE X2AVE[i]
            x2ave[i] = x2ave[i] + (currentPosition * currentPosition);
            cout << "==== The squared distance is <x^2> = " << x2ave[i] << endl << endl;
        }
        cout << endl;
    }

    //Mean Square
    cout << endl;
    for (int i = 0; i < numSteps; i++) {
        double meanSquare = x2ave[i] / numWalkers;
        posMean << i << ", " << meanSquare << endl;
    }

    //PROBABILITY
    cout << endl;
    cout << "Probability Values:" << endl << endl;
    for (int i = 0; i < (2 * numSteps + 1); i++)
    {
        //COMPUTE PROBABILITY OF POSITION i
        prb[i] = (double)countNumberOfVisits[i] / (double)numWalkers;
        if (countNumberOfVisits[i] != 0) {
            cout << "Total Steps and Probability of Walker To Be At Distance " << i - numSteps << ": " << countNumberOfVisits[i] << ", " << prb[i] << endl;
        }
        probOut << i - numSteps << ", " << prb[i] << endl;  // Output(x,y): POSITION VS PROBABILITY
        posNSteps << i - numSteps << ", " << countNumberOfVisits[i] << endl; // Output(x,y): POSITION VS TOTAL STEPS THE POSITION WAS VISITED
    }

    //THIS IS JUST A REMINDER OF THE THREE OUTPUT FILES CREATED
    cout << endl
        << "3 Output files have been created: " << endl
        << "positionVSProbabilityOutput.csv, " << endl
        << "positionVSNumStepsOutput.csv, " << endl
        << "positionVSMeanSquareOutput.csv," << endl << endl;

    probOut.close();
    posNSteps.close();
    posMean.close();
    return 0;
}
