#pragma once

#include <iostream>
#include <fstream>
#include <string>

#include "parameters.hpp"

using namespace std;

class Reader
{
public:
	Reader() {};
	~Reader() {};
	pair<double, double> * readData(size_t i); //returns pointer to data array
	void readDataLoc(); // to read data location
	void readResLoc(); // to read where to put the results
	void readParams(); // to read parameters
	string getResPath(); // to get location path
	void getParams(size_t &minChirpSize, size_t &minCandSize, size_t &ptsForCorrelCheck, double &stepSizeForCand, double &maxDistBetweenPts,
     double &EPS_start, double &EPS_end, double &maxDistBetweenCands, double &EPS_angle, double &R0); // to get parameters
	size_t getLength(size_t i); // returns length of ith data
	size_t nProc(); // returns total number of files to be rpocessed
protected:
	string _dataLocationT[__MAX_FILES_TO_PROCESS__];
	string _dataLocationF[__MAX_FILES_TO_PROCESS__];
	string _resLocation;
	size_t _dataLength[__MAX_FILES_TO_PROCESS__];
	size_t _nFilesToProcess;
	size_t _minChirpSize;
	size_t _minCandSize;
	size_t _ptsForCorrelCheck;
	double _stepSizeForCand;
	double _maxDistBetweenPts;
	double _EPS_start;
	double _EPS_end;
	double _maxDistBetweenCands;
	double _EPS_angle;
	double _R0;
};

pair<double, double> * Reader::readData(size_t i)
{
	pair<double, double> * data = new pair<double, double>[_dataLength[i]];

	ifstream dat_t, dat_f;

	dat_f.open(_dataLocationF[i]);
	if (!dat_f)
	{
		cout << "Can't open " << _dataLocationF[i] << endl;
	}
	else
	{
		for (int j = 0; j < _dataLength[i]; ++j)
		{
			dat_f >> data[j].first;
		}
		cout << _dataLength[i] << " elements of " << _dataLocationF[i] << " are read" << endl;
	}
	dat_f.close();
	dat_t.open(_dataLocationT[i]);
	if (!dat_t)
	{
		cout << "Can't open " << _dataLocationT[i] << endl;
	}
	else
	{
		for (int j = 0; j < _dataLength[i]; ++j)
		{
			dat_t >> data[j].second;
		}
		cout << _dataLength[i] << " elements of " << _dataLocationT[i] << " are read" << endl;
	}
	
	dat_t.close();
	return data;
}

void Reader::readParams()
{
	ifstream par;
	par.open("parameters.txt");
	if (!par)
	{
		cout << endl << "Can't open parameters.txt, default parameters will be set" << endl;
		_minChirpSize = __MIN_CHIRP_SIZE__;
		_minCandSize = __MIN_CAND_SIZE__;
		_ptsForCorrelCheck = __PTS_FOR_CORREL_CHECK__;
		_stepSizeForCand = __STEP_SIZE_FOR_CAND__;
		_maxDistBetweenPts = __MAX_DIST_BETWEEN_POINTS__;
		_EPS_start = __EPS_START__;
		_EPS_end = __EPS_END__;
		_maxDistBetweenCands = __MAX_DIST_BETWEEN_CANDS__;
		_EPS_angle = __EPS_ANGLE__;
		_R0 = __R0__;
	}
	else
	{
		par >> _minChirpSize;
		par >> _minCandSize;
		par >> _ptsForCorrelCheck;
		par >> _stepSizeForCand;
		par >> _maxDistBetweenPts;
		par >> _EPS_start;
		par >> _EPS_end;
		par >> _maxDistBetweenCands;
		par >> _EPS_angle;
		par >> _R0;
	}
	par.close();
	cout << endl << "Algorithm will start with the following parameters:" << endl << endl << " minChirpSize = " << _minChirpSize << endl << " minCandSize = "
	 << _minCandSize << endl << " ptsForCorrelCheck = " << _ptsForCorrelCheck << endl << " stepSizeForCand = " << _stepSizeForCand << endl <<
      " maxDistBetweenPts = " << _maxDistBetweenPts << endl << " EPS_start = " << _EPS_start <<  endl << " EPS_end = " << _EPS_end << endl << " maxDistBetweenCands = "
	  << _maxDistBetweenCands  <<  endl << " EPS_angle = " << _EPS_angle << endl << " R0 = " << _R0 << endl << endl;
}

void Reader::readResLoc()
{
	ifstream res_loc;

	res_loc.open("res_loc.txt");
	if (!res_loc)
	{
		cout << "Can't open res_loc.txt" << endl;
	}
	else
	{
		res_loc >> _resLocation;
	}
	res_loc.close();
}

void Reader::readDataLoc()
{
	ifstream dat_loc;
	size_t n = 0;

	dat_loc.open("data_loc.txt");
	if (!dat_loc)
	{
		cout << "Can't open data_loc.txt" << endl;
	}
	else
	{
		while (dat_loc >> _dataLocationT[n])
		{
			dat_loc >> _dataLocationF[n];
			dat_loc >> _dataLength[n];
			++n;
		}
	}
	dat_loc.close();
	_nFilesToProcess = n;
}

void Reader::getParams(size_t &minChirpSize, size_t &minCandSize, size_t &ptsForCorrelCheck, double &stepSizeForCand, double &maxDistBetweenPts,
     double &EPS_start, double &EPS_end, double &maxDistBetweenCands, double &EPS_angle, double &R0)
{
	minChirpSize = _minChirpSize;
	minCandSize = _minCandSize;
	ptsForCorrelCheck = _ptsForCorrelCheck;
	stepSizeForCand = _stepSizeForCand;
	maxDistBetweenPts = _maxDistBetweenPts;
	EPS_start = _EPS_start;
	EPS_end = _EPS_end;
	maxDistBetweenCands = _maxDistBetweenCands;
	EPS_angle = _EPS_angle;
	R0 = _R0;	
}

size_t Reader::nProc()
{
	return _nFilesToProcess;
}

size_t Reader::getLength(size_t i)
{
	return _dataLength[i];
}

string Reader::getResPath()
{
	return _resLocation;
}