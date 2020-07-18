#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <bits/stdc++.h> 

using namespace std;

const int numWalkers = 1000;  //NUMBER OF WALKERS
const int numSteps = 22;  //NUMBER OF STEPS (CANNOT SET GREATER THAN 22 SINCE FACTORIAL IS GOING TO BEHAVE ABNORMALLY, MEMORY CAN'T CONTAIN)

double q1 = 1 / sqrt(numSteps), q2 = 1, q3 = sqrt(numSteps);

unsigned long long fact(int n)
{
    if ((n == 0) || (n == 1))
        return 1;
    else
        return n * fact(n - 1);
}

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
    double* cosRq1 = new double[numSteps];
    double* cosRq2 = new double[numSteps];
    double* cosRq3 = new double[numSteps];

    for (int itime = 0; itime < numSteps; itime++) {
        aStepSquared[itime] = 0;
        cosRq1[itime] = 0.0;
        cosRq2[itime] = 0.0;
        cosRq3[itime] = 0.0;
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
            //TASK TO SHOW NUMERICAL AND ANALYTICAL PART OF x^2 or aStepSquared
            cout << "==== New position is: " << currentPosition << endl;

            //COMPUTE aStepSquared[itime] and cosRq1, cosRq2, cosRq3
            //TASK TO SHOW BOTH NUMERICAL AND ANALYTIC PART OF THE DIFFERENT q's: q1, q2, q3
            aStepSquared[itime] += (currentPosition * currentPosition);
            cosRq1[itime] += cos(q1 * currentPosition);
            cosRq2[itime] += cos(q2 * currentPosition);
            cosRq3[itime] += cos(q3 * currentPosition);

            //cout << "==== The squared displacement is x^2 = " << aStepSquared[itime] << endl << endl;

        }
        cout << endl;
        //INCREMENT countNumberOfVisits[] EVERYTIME currentPosition IS VISITED
        countNumberOfVisits[currentPosition + numSteps] += 1;
    }

    //Mean Square
    cout << endl;
    for (int itime = 0; itime < numSteps; itime++) {
        double meanSquare = (double)aStepSquared[itime] / (double)numWalkers;
        cosRq1[itime] = (double)cosRq1[itime] / (double)numWalkers;
        cosRq2[itime] = (double)cosRq2[itime] / (double)numWalkers;
        cosRq3[itime] = (double)cosRq3[itime] / (double)numWalkers;

        posMean << itime << ", " << meanSquare << ", " << cosRq1[itime] << ", " << cosRq2[itime] << ", " << cosRq3[itime] << endl;
    }


    //GETTING PROBABILITY OF FINAL POSITION
    //USING FACTORIAL METHOD; ANALYTICAL VERSION
    //TASK: TO INSERT AGAIN THE PREVIOUS PROBABILITY AND OVERLAP THE VALUES
    //TO SHOW PLOTS OF PROBABILITY DIST: GAUSSIAN, BINOMIAL (USE EXCEL), NUMERICAL SOLUTION
    cout << endl;
    cout << "Probability Values:" << endl << endl;

    for (int position = (-1 * numSteps); position <= numSteps; position+=2)
    {
                double totalPossibleVisits = fact(numSteps) / (fact((numSteps + abs(position)) / 2) * fact((numSteps - abs(position)) / 2));
                double prob = totalPossibleVisits / pow(2, numSteps);
                probOut << position << ", " << prob << endl;
    }

    cout << endl << endl;

    probOut.close();
    posMean.close();
    return 0;
}