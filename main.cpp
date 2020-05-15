#include <iostream>

#include "diplo.hpp"
#include "reader.hpp"

using namespace std;

int main()
{
	Reader reader; //needed to read parameters and res/data location
	Algorithm alg; // algorithm implementation and writing results
	
	/* Parameters */
	size_t nProc;
	size_t curLength;
	size_t minChirpSize;
	size_t minCandSize;
	size_t ptsForCorrelCheck;
	double stepSizeForCand;
	double maxDistBetweenPts;
	double EPS_start;
	double EPS_end;
	double maxDistBetweenCands;
	double EPS_angle;
	double R0;

	double time;

	/* Reading parameters, data and result location */
	reader.readParams(); //same parameters for all data sets 
	reader.readDataLoc();
	reader.readResLoc();
	reader.getParams(minChirpSize, minCandSize, ptsForCorrelCheck, stepSizeForCand, maxDistBetweenPts, EPS_start, EPS_end, maxDistBetweenCands, EPS_angle, R0);
	nProc = reader.nProc();
 //possible to add parallelezation
	for (size_t i = 0; i < nProc; ++i){

		cout << endl << "Processing " << i << " data" << endl;

		curLength = reader.getLength(i);

		/* Processing data */
		/* Setting input data for the algorithm */
		alg.init(reader.readData(i), curLength, minChirpSize, minCandSize, ptsForCorrelCheck, stepSizeForCand, maxDistBetweenPts, EPS_start, EPS_end,
		 maxDistBetweenCands, EPS_angle, R0);
		
		/* Starting the algorithm */
		alg.run();

		/* Writing results and statistics */
		alg.writeStats(i, reader.getResPath());
		alg.writeResults(i, reader.getResPath());

		time = alg.getTime(); // in seconds
		cout << "Execution time for " << i << " data: " << time << " seconds" << endl;

		/* Removing initialized data */
		alg.clear();
	}

	system("pause");
	return 0;
}