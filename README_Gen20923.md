#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <chrono>
#include <random>

using namespace std;
// Use .dat instead of .txt

int main(int argc, char const* argv[]) {

	/* generate numWalkers and numSteps between 1 and max*/
	int numWalkers = 1;
	int numSteps = 20;

	float* stepsTracker = new float[numSteps * numWalkers]; //tracks all steps taken by walkers. 
	for (int i = 0; i < numSteps * numWalkers; i++)
		stepsTracker[i] = 0; //initialize tracker

	ofstream position; //output stream for all final steps taken by walkers
	//http://www.cplusplus.com/doc/tutorial/files/
	fstream positionBinary("positions2.bin", ios::out | ios::app | ios::binary); //output stream for binary
	position.open("position2.txt"); //TASK to use position instead of the long name

	// the following FOR LOOP prints all steps taken by each walker
	// OUTPUT FORMAT IS  :
	// walker1 wlkrStep1 wlkrStep2 ...
	// walker2 wlkrStep1 wlkrStep2 ...
	// walker3 wlkrStep1 wlkrStep2 ...

	for (int walker = 0; walker < numWalkers; walker++) {

		double currentPosition = 0.0; //Change depending on the current position of walker

		for (int step = 0; step < numSteps; step++) {

			//HEAD OR TAIL CONCEPT
			//https://www.cplusplus.com/reference/random/normal_distribution/normal_distirbution/

			//Gaussian or normal random number
			unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
			std::default_random_engine generator(seed);
			std::normal_distribution<double> distribution(0.0, 1.0);

			double randomNum = distribution(generator);
			cout << "randomNum " << randomNum ;

			currentPosition += randomNum;
			cout << "  currentPosition " << currentPosition << endl;

			stepsTracker[walker * numSteps + step] = currentPosition;
		}
	}

	// this FOR LOOP prints all final steps of each walker per time or we use "step" here
	// we will utilize the stepsTracker values in the previous loop to fetch all final steps
	// OUTPUT FORMAT IS  :
	// time1 wlkr1 wlkr2 wlkr3 ...
	// time2 wlkr1 wlkr2 wlkr3 ...
	// time3 wlkr1 wlkr2 wlkr3 ...

	for (int step = 0; step < numSteps; step++) {

		int theStep = step + 1;
		position << theStep;
		positionBinary.write((char*)&theStep, sizeof(int));

		//FOR LOOP FOR EACH WALKER
		for (int walker = 0; walker < numWalkers; walker++) {

			//TASK use wlkr * numSteps + step instead of trackerCounter
			double finalStep = stepsTracker[walker * numSteps + step];
			position << " " << finalStep;

			positionBinary.write((char*)&finalStep, sizeof(int));
		}
		position << endl;
	}

	position.close();
	positionBinary.close();

	return 0;
}
