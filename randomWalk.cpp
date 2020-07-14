#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <string>

using namespace std;

//CONSTANT VARIABLES
const int numberOfWalkers = 2;
const int numSteps = 7;
int accum = 0;

//START MAIN PROGRAM
int main(int argc, char const* argv[])
{
    cout << "NUMBER OF WALKERS: " << numberOfWalkers << endl;
    cout << "NUMBER OF STEPS: " << numSteps << endl << endl;

    //INITIALIZING FOR TEXT FILES
    ofstream probOut;
    ofstream posNSteps;
    ofstream posMean;
    ofstream stepsX2ave;
    ofstream probCounter;

    probOut.open("positionVSProbabilityOutput.csv");
    posNSteps.open("positionVSNumStepsOutput.csv");
    posMean.open("positionVSMeanSquareOutput.csv");
    stepsX2ave.open("x2averageVSsteps.csv");
    probCounter.open("myCounterVSprobability.csv");

    //DECLARATION OF VARIABLES & ARRAYS

      // x2ave IS THE MEAN SQUARE DISPLACEMENT AND IS A FUNCTION OF NUMBER OF STEPS OR TIME, VARIABLE i
      // x IS DISPLACEMENT AFTER A SINGLE WALKER IS DONE AND IS A FUNCTION OF NUMBER OF WALKERS, VARIABLE k
    int* step_number = new int[numSteps];
    int* x2ave = new int[numSteps];
    int* x = new int[numberOfWalkers];


    //  TO BE USED IN SOLVING FOR PROBABILITY DISTRIBUTION
    double* prb = new double[2 * numSteps + 1];
    int* myCounter = new int[numSteps * 2 + 1];

    // DECLARATION OF ARRAY AS ZEROS
    for (int i = 0; i < numSteps; i++) {
        x2ave[i] = 0;
        x[i] = 0;
        step_number[i] = 0;
    }

    for (int i = 0; i < (2 * numSteps + 1); i++) {
        myCounter[i] = 0;
        prb[i] = 0.0;
    }

    // LOOP FOR A NUMBER OF WALKERS TAKING SUCH NUMBER OF RANDOM STEPS
    for (int k = 0; k < numberOfWalkers; k++) {

        int position = 0; //SET INITIAL POSITION TO ZERO PER WALKER K = 1,2,3,.., numberOfWalkers

        for (int i = 0; i < numSteps; i++)  // START OF LOOP FOR EVERY RANDOM STEP OF A SINGLE WALKER
        {
            int randomNum = rand() % 2;

            if (randomNum == 1)
            {
                position++;
            }
            else
            {
                position--;
            }

            //now you get a new position 
            x2ave[i] = x2ave[i] + position * position;
            cout << "x2ave[i=" << i << "] " << x2ave[i] << endl << endl;
            accum = accum + x2ave[i]; //ADDS UP ALL THE X2AVE

        }
        cout << "New position for a single walker:  " << position << endl;
        x[k] = position;    // x[k] IS THE POSITION WHERE THE WALKER k STOPPED WALKING
        cout << "x[" << k << "] at k = " << k << " with upper limit " << numberOfWalkers << " is ===> " << x[k] << endl;
        myCounter[x[k] + numSteps] += 1;    // TO BE USED IN PLOTTING THE PROBABILITY DISTRIBUTION
        cout << "myCounter[" << x[k] + numSteps << "] is ===> " << myCounter[x[k] + numSteps] << endl;
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++ " << endl << endl;

        // NOW WE START THE LOOP FOR VARIABLE k FOR A NEW WALKER UNTIL WE'RE DONE FOR ALL WALKERS

    }// END OF THE LOOP FOR A NUMBER OF WALKERS TAKING SUCH NUMBER OF RANDOM STEPS


    //MEAN SQUARE
    cout << endl;
    for (int i = 0; i < numSteps; i++) {
        cout << "x2ave[i] " << x2ave[i] << endl;
        double meanSquare = x2ave[i] / numberOfWalkers;
        cout << i << "," << meanSquare << endl;
        posMean << i << "," << meanSquare << endl; // Output: meanSquare vs numsteps
    }

    probCounter.close();
    probOut.close();
    posNSteps.close();
    posMean.close();
    stepsX2ave.close();
    return 0;
}
