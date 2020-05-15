#pragma once

#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <stack>
#include <iomanip>

#include "external.hpp"

using namespace std;

class Algorithm
{
public:
	Algorithm(){};
	~Algorithm(){};
	void init(pair<double, double>* data, size_t dataLength, size_t &minChirpSize, size_t &minCandSize, size_t &ptsForCorrelCheck, double &stepSizeForCand,
	 double &maxDistBetweenPts, double &EPS_start, double &EPS_end, double &maxDistBetweenCands, double &EPS_angle, double &R0); // to initialize algorithm parameters and data
	void run(); // to start algorithm
	void clear(); // to reset parameters
	double getTime(); // to get execution time in seconds
	void writeResults(size_t k, const string res = ""); // to write results
	void writeStats(size_t i, const string sts = ""); // to write current parameters
protected:
	void process();
	void prolongMethod (size_t start);
	void mergeCandidates(size_t start);
	vector<pair<double, double>> add(vector<pair<double, double>> first, vector<pair<double, double>> second);
protected:
	size_t _dataLength;
	size_t _mx;
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
	double _execTime;
	pair<double, double>* _data;
	vector<vector<pair<double, double>>> _chirpCandidates;
	vector<vector<pair<double, double>>> _chirps;
	vector <double> _cand_time_beg;
	vector <double> _cand_time_end;
	bool* _marks;
	bool* _marksCand;
};

void Algorithm::init(pair<double, double>* data, size_t dataLength, size_t &minChirpSize, size_t &minCandSize, size_t &ptsForCorrelCheck, double &stepSizeForCand,
	 double &maxDistBetweenPts, double &EPS_start, double &EPS_end, double &maxDistBetweenCands, double &EPS_angle, double &R0)
{
	_data = data;
	_marks = new bool[dataLength];
	for (size_t i = 0; i < dataLength; ++i)
	{
		_marks[i] = false;
	}
	_dataLength = dataLength;
	_minChirpSize = minChirpSize;
	_minCandSize = minCandSize;
	_ptsForCorrelCheck = ptsForCorrelCheck;
	_stepSizeForCand = stepSizeForCand;
	_maxDistBetweenPts = maxDistBetweenPts;
	_EPS_start = EPS_start;
	_EPS_end = EPS_end;
	_maxDistBetweenCands = maxDistBetweenCands;
	_EPS_angle = EPS_angle;
	_R0 = R0;
	_mx = 0;
}

void Algorithm::run()
{
	double begin = (double) clock();
	sort(_data, _data + _dataLength, pairComparison);
	process();
	double end = (double) clock();

	cout << _dataLength << " elements processed." << endl << endl;

	_execTime = (double)(end - begin) / (double) CLOCKS_PER_SEC;
}

void Algorithm::process()
{
	/* Searching candidates */
	for (size_t i = _dataLength - 1; i > 0; --i)
	{
		if (!_marks[i]) {
			prolongMethod(i);
		}
	}
	delete _data;
	delete[] _marks;
	cout << "Candidates size = " << _chirpCandidates.size() << endl;

	_marksCand = new bool[_chirpCandidates.size()];
	for (size_t i = 0; i < _chirpCandidates.size(); ++i)
	{
		_marksCand[i] = false;
	}

	/* Merging candidates */
	sort(_chirpCandidates.begin(), _chirpCandidates.end(), vectorComparison);
	for (size_t i = 0; i < _chirpCandidates.size(); ++i)
	{
		_cand_time_beg.push_back(_chirpCandidates[i][0].second);
		_cand_time_end.push_back(_chirpCandidates[i][_chirpCandidates[i].size() - 1].second);
	}
	for (size_t i = 0; i < _chirpCandidates.size(); ++i)
	{
		if (!_marksCand[i])
		{
			mergeCandidates(i);
		}
	}
	sort(_chirps.begin(), _chirps.end(), vectorComparison);
	cout << "Chirps size = " << _chirps.size() << endl;
	cout << "The longest chirp contains " << _mx << " elements" << endl;
	delete [] _marksCand;
}

void Algorithm::clear()
{
	_dataLength = 0;
	_minChirpSize = 0;
	_minCandSize = 0;
	_ptsForCorrelCheck = 0;
	_stepSizeForCand = 0.0;
	_maxDistBetweenPts = 0.0;
	_EPS_start = 0.0;
	_EPS_end = 0.0;
	_maxDistBetweenCands = 0.0;
	_EPS_angle = 0.0;
	_R0 = 0.0;
	_execTime = 0.0;
	_chirps.clear();
	_cand_time_end.clear();
	_cand_time_beg.clear();
	_chirpCandidates.clear();
}

void Algorithm::writeStats(size_t i, string sts)
{
	string num = std::to_string(i);
	ofstream stats;
	stats.open(sts + "\\stats_" + num + ".txt");

	stats << "Data was processed with the following parameters: " << endl << endl;
	stats << "minChirpSize: " << _minChirpSize << endl;
	stats << "minCandSize: " << _minCandSize << endl;
	stats << "ptsForCorrelCheck: " << _ptsForCorrelCheck << endl;
	stats << "stepSizeForCand: " << _stepSizeForCand << endl;
	stats << "maxDistBetweenPts: " << _maxDistBetweenPts << endl;
	stats << "EPS_start: " << _EPS_start << endl;
	stats << "EPS_end: " << _EPS_end << endl;
	stats << "maxDistBetweenCands: " << _maxDistBetweenCands << endl;
	stats << "EPS_angle: " << _EPS_angle << endl;
	stats << "R0: " << _R0 << endl << endl;

	stats << "Number of chirps observed: " << _chirps.size() << endl;
	stats << "The longest chirp observed contains " << _mx << " elements" << endl;

}

void Algorithm::writeResults(size_t k, const string res)
{
	string num = std::to_string(k);
	ofstream dat;
	dat.open(res + "chirps_f_" + num + ".txt");
	dat.precision(10);
	for (size_t i = 0; i < _chirpCandidates.size(); ++i)
	{
		for (size_t j = 0; j < _chirpCandidates[i].size(); ++j)
		{
			dat << _chirpCandidates[i][j].first << ";";
		}
		dat << "\n";
	}
	dat.close();
	dat.open(res + "chirps_t_" + num + ".txt");
	dat.precision(10);
	for (size_t i = 0; i < _chirpCandidates.size(); ++i) 
	{
		for (size_t j = 0; j < _chirpCandidates[i].size(); ++j)
		{
			dat << _chirpCandidates[i][j].second << ";";
		}
		dat << "\n";
	}
	dat.close();
}

double Algorithm::getTime()
{
	return _execTime;
}

void Algorithm::prolongMethod (size_t start)
{
	vector<pair<double, double>> temp;
	pair<double, double> currentPair;
	vector <double> localDer;
	temp.push_back(_data[start]);
	int s = 2;
	size_t n_prolong = (size_t) (_maxDistBetweenPts / _stepSizeForCand);
	int numOfAddedPts = 0;

	/* Adding first _ptsForCorrelCheck points for further MSE */
	while (s > 0 && temp.size() <= _ptsForCorrelCheck)
	{
		--s;
		if ((s % 2) == 1)
		{
			currentPair = temp[temp.size() - 1];
		}
		else{
			currentPair = temp[0];
		}

		pair<double, double> temppUp = {_maxDistBetweenPts, _maxDistBetweenPts};
		pair<double, double> temppDn = {- _maxDistBetweenPts, - _maxDistBetweenPts};
		temppUp.first += currentPair.first;
		temppDn.first += currentPair.first;
		temppUp.second += currentPair.second;
		temppDn.second += currentPair.second;
		size_t up = upper_bound(&_data[0], &_data[_dataLength - 1], temppUp, pairComparison) - &_data[0];
		size_t dn = lower_bound(&_data[0], &_data[_dataLength - 1], temppDn, pairComparison) - &_data[0];
		if (dn > 0) --dn;
		if (up < (_dataLength - 1)) ++up;
		/* Adding new points next to the current one */
		for (size_t i = up; i > dn; --i)
		{
			if ((_marks[i] == false) && (distance(currentPair, _data[i]) < (_maxDistBetweenPts))) 
			{
				_marks[i] = true;
				temp.push_back(_data[i]);
				(++s)++;
			}
		}
		while ((temp.size() > _ptsForCorrelCheck) && (!isCorrelated(&temp, 0, temp.size() - 1, _R0)))
		{
			temp.pop_back();
			--s;
		}
		sort(temp.begin(), temp.end(), pairComparison);

	}

	if (temp.size() > _ptsForCorrelCheck && isCorrelated(&temp, 0, temp.size() - 1, _R0))
	{
		int check = 1;
		pair <double, double> coeffs;
		vector <size_t> localMarks;
		int ctr = 0;
		pair <double, double> virtualPoint;
		size_t l =1;
		double d = 0;
		size_t corrCheck = _ptsForCorrelCheck;
		/* Adding points from back */
		while (check > 0 && (l < _dataLength))
		{
			++l;
			--check;
			coeffs = linearRegr(&temp, temp.size() - 1 - _ptsForCorrelCheck, temp.size() - 1);
			for (size_t i = 1; i <= n_prolong; ++i)
			{
				virtualPoint = {coeffs.first * (temp[temp.size() - 1].second + i * _stepSizeForCand) + coeffs.second, temp[temp.size() - 1].second + i * _stepSizeForCand};
				double coeff = _stepSizeForCand * (_EPS_start + sqrt((double) i / ((double) (n_prolong))) * (_EPS_end - _EPS_start)); //   _maxDistBetweenPts * (0.25 + sqrt((double) i / ((double) (n_prolong) * 0.8)));
				pair<double, double> temppUp = {coeff, coeff};
				pair<double, double> temppDn = {-coeff, -coeff};
				temppDn.first += virtualPoint.first;
				temppUp.first += virtualPoint.first;
				temppDn.second += virtualPoint.second;
				temppUp.second += virtualPoint.second;
				size_t up = upper_bound(&_data[0], &_data[_dataLength - 1], temppUp, pairComparison) - &_data[0];
				size_t dn = lower_bound(&_data[0], &_data[_dataLength - 1], temppDn, pairComparison) - &_data[0];
				if (up < (_dataLength - 1)) ++up;
				if (dn > 0) --dn;
				/* Adding new points next to the current one */
				for (size_t j = dn; j < up; ++j)
				{
					if (!_marks[j] &&(distance(virtualPoint, _data[j]) < coeff))// && (temp[temp.size() - 1].second < _data[j].second))
					{
						_marks[j] = true;
						temp.push_back(_data[j]);
						localMarks.push_back(i);
						++check;
						++numOfAddedPts;
						break;
					}
				}
				if (check > 0) break;
			}
		}
		//if (numOfAddedPts > (_ptsForCorrelCheck - 3))
		//{
		//	numOfAddedPts = _ptsForCorrelCheck - 3;
		//}
		//for (int i = temp.size() - 1; i > temp.size() - numOfAddedPts; --i)
		//{
		//	double x = temp[i].second - temp[i - 1].second;
		//	double y = temp[i].first - temp[i - 1].first;
		//	localDer.push_back(y / x);
		//}
		//for (int i = localDer.size() - 1; i > 1; --i)
		//{
		//	if (fabs(localDer[i] - localDer[i - 1] > 0.6))
		//	{
		//		ctr++;
		//	}
		//}
		//if (ctr > 1)
		//{
		//	for (int i = localMarks.size() - 1; i > localMarks.size() - ctr; --i)
		//	{
		//		_marks[localMarks[i]] = false;
		//		temp.pop_back();
		//	}
		//}
		//ctr = 0;
		//for (int i = temp.size() - 1; i > temp.size() - numOfAddedPts; --i)
		//{
		//	if(distance(temp[i], temp[i - 1]) > (_maxDistBetweenPts))
		//	{
		//		ctr++;
		//	}
		//}
		//for (int i = temp.size() - 1; i > temp.size() - numOfAddedPts; --i)
		//{
		//	if(distance(temp[i], temp[i - 1]) > (_maxDistBetweenPts))
		//	{
		//		ctr++;
		//	}
		//}
		//if (ctr > 2)
		//{
		//	cout << localMarks.size() << endl;
		//	for (size_t i = localMarks.size(); i > localMarks.size() - ctr; --i)
		//	{
		//		cout << localMarks[i - 1] << endl;
		//		_marks[localMarks[i - 1]] = false;
		//		temp.pop_back();
		//	}
		//}
		//if (corrCheck > localMarks.size()) corrCheck = localMarks.size();
		//if (!isCorrelated(&temp, temp.size() - corrCheck, temp.size() - 1, _R0))
		//{
		//	for (size_t k = 0; k < corrCheck; ++k)
		//	{
		//		_marks[localMarks[localMarks.size() - k - 1]] = false;
		//		temp.pop_back();
		//	}
		//}

		corrCheck = _ptsForCorrelCheck;
		ctr = 0;
		numOfAddedPts = 0;
		localMarks.clear();
		localDer.clear();
		check = 1;
		l = 1;
		d = 0;
		/* Adding points from head */
		while (check > 0 && (l < _dataLength))
		{
			++l;
			--check;
			coeffs = linearRegr(&temp, 0, _ptsForCorrelCheck - 1);
			for (size_t i = 1; i <= n_prolong; ++i)
			{
				virtualPoint = {coeffs.first * (temp[0].second - i * _stepSizeForCand) + coeffs.second, temp[0].second - i * _stepSizeForCand};
				//cout << "From head: virtualPoint = {" << virtualPoint.second << ", " << virtualPoint.second << "} i = " << i <<endl;
				double coeff = _stepSizeForCand * _EPS_start + sqrt((double) i / ((double) (n_prolong))) * (_EPS_end - _EPS_start) * _stepSizeForCand; //_maxDistBetweenPts * (0.25 + sqrt((double) i / ((double) (n_prolong) * 0.8)));
				pair<double, double> temppUp = {coeff, coeff};
				pair<double, double> temppDn = {-coeff, -coeff};
				temppDn.first += virtualPoint.first;
				temppUp.first += virtualPoint.first;
				temppDn.second += virtualPoint.second;
				temppUp.second += virtualPoint.second;
				size_t up = upper_bound(&_data[0], &_data[_dataLength - 1], temppUp, pairComparison) - &_data[0];
				size_t dn = lower_bound(&_data[0], &_data[_dataLength - 1], temppDn, pairComparison) - &_data[0];
				if (up < (_dataLength - 1)) ++up;
				if (dn > 0) --dn;
				/* Adding new points next to the current one */
				for (size_t j = up; j > dn; --j)
				{
					if (!_marks[j] && (distance(virtualPoint, _data[j]) < coeff))// && (temp[0].second > _data[j].second))
					{
						_marks[j] = true;
						temp.push_back(_data[j]);
						localMarks.push_back(i);
						++check;
						++numOfAddedPts;
						break;
					}
				}
				sort(temp.begin(), temp.end(), pairComparison);
				if (check > 0) break;
			}
		}
		//if (numOfAddedPts > (_ptsForCorrelCheck - 3))
		//{
		//	numOfAddedPts = _ptsForCorrelCheck - 3;
		//}
		vector <pair <double, double>> revTmp;
		for (size_t i = 1; i <= temp.size(); ++i)
		{
			revTmp.push_back(temp[temp.size() - i]);
		}
		//for (size_t i = 1; i < numOfAddedPts; ++i)
		//{
		//	double x = temp[i].second - temp[i-1].second;
		//	double y = temp[i].first - temp[i-1].first;
		//	localDer.push_back(y / x);
		//}
		//for (int i = localDer.size() - 1; i > 1; --i)
		//{
		//	if (fabs(localDer[i] - localDer[i - 1] > 0.6))
		//	{
		//		ctr++;
		//	}
		//}
		//if (ctr > 1)
		//{
		//	for (int i = localMarks.size() - 1; i > localMarks.size() - ctr; --i)
		//	{
		//		_marks[localMarks[i]] = false;
		//		revTmp.pop_back();
		//	}
		//}
		//ctr = 0;
		//for (size_t i = 1; i <= numOfAddedPts; ++i)
		//{
		//	if(distance(temp[i], temp[i - 1]) > (_maxDistBetweenPts))
		//	{
		//		ctr++;
		//	}
		//}
		//if (ctr > 2)
		//{
		//	for (size_t i = localMarks.size() - 1; i > localMarks.size() - ctr; --i)
		//	{
		//		_marks[localMarks[i]] = false;
		//		revTmp.pop_back();
		//	}
		//}
		//if (corrCheck > localMarks.size()) corrCheck = localMarks.size();
		//if (!isCorrelated(&revTmp, revTmp.size() - corrCheck, revTmp.size() - 1, _R0))
		//{
		//	for (size_t k = 0; k < corrCheck; ++k)
		//	{
		//		_marks[localMarks[localMarks.size() - k - 1]] = false;
		//		revTmp.pop_back();
		//	}
		//}
		if (revTmp.size() > _ptsForCorrelCheck + 1)
		{
			sort(revTmp.begin(), revTmp.end(), pairComparison);
			_chirpCandidates.push_back(revTmp);
		}
	}
	temp.clear();

}

void Algorithm::mergeCandidates(size_t start)
{
	size_t st = start;
	if (!_marksCand[start])
	{
		_marksCand[start] = true;
		vector<pair<double, double>> tmpVec = _chirpCandidates[start];
		bool possible_to_add = true;
		while(possible_to_add)
		{
			possible_to_add = false;

			/* Trying to add from back */
			double temppup = _maxDistBetweenCands;
			double temppdn = - _maxDistBetweenCands;
			temppup += (tmpVec[tmpVec.size() - 1]).second;
			temppdn += (tmpVec[tmpVec.size() - 1]).second;
			size_t up = upper_bound(_cand_time_beg.begin(), _cand_time_beg.end(), temppup) - _cand_time_beg.begin();
			if (up < (_chirpCandidates.size() - 1)) ++up;
		
			for (size_t j = st; j < up; ++j)
			{
				pair <double, double> beg_der = linearRegr(&_chirpCandidates[j], 0, _ptsForCorrelCheck + 1);
				pair <double, double> end_der = linearRegr(&tmpVec, tmpVec.size() - _ptsForCorrelCheck - 2, tmpVec.size() - 1);
				pair <double, double> predicted_pt = {_chirpCandidates[j][0].second * end_der.first + end_der.second, _chirpCandidates[j][0].second};
				if ((fabs(1 - end_der.first / beg_der.first) < _EPS_angle) && (distance(predicted_pt, _chirpCandidates[j][0]) < (_maxDistBetweenPts))
				 && (tmpVec[tmpVec.size() - 1].second < _chirpCandidates[j][0].second) && !_marksCand[j])
				{
					tmpVec = add(tmpVec, _chirpCandidates[j]);
					_marksCand[j] = true;
					possible_to_add = true;
				}
			}
			
			/* Trying to add from front */
			temppup = _maxDistBetweenCands;
			temppdn = - _maxDistBetweenCands;
			temppup += (tmpVec[0]).second;
			temppdn += (tmpVec[0]).second;
			size_t loc_st = st;
			size_t dn = lower_bound(_cand_time_end.begin(), _cand_time_end.end(), temppdn) - _cand_time_end.begin();		
			for (size_t j = dn; j < st; ++j)
			{
				pair <double, double> beg_der = linearRegr(&tmpVec, 0, _ptsForCorrelCheck + 1);
				pair <double, double> end_der = linearRegr(&_chirpCandidates[j], _chirpCandidates[j].size() - _ptsForCorrelCheck - 2, _chirpCandidates[j].size() - 1);
				pair <double, double> predicted_pt = {_chirpCandidates[j][_chirpCandidates[j].size() - 1].second * beg_der.first + beg_der.second, _chirpCandidates[j][_chirpCandidates[j].size() - 1].second};
				if ((fabs(1 - end_der.first / beg_der.first) < _EPS_angle) && (distance(predicted_pt, _chirpCandidates[j][_chirpCandidates[j].size() - 1]) < (_maxDistBetweenPts))
				 && (tmpVec[0].second > _chirpCandidates[j][_chirpCandidates[j].size() - 1].second) && !_marksCand[j])
				{
					tmpVec = add(_chirpCandidates[j], tmpVec);
					_marksCand[j] = true;
					possible_to_add = true;
					loc_st = j;
				}
			}
			st = loc_st;
		}
		if (tmpVec.size() > _minChirpSize)
		{
			if (_mx < tmpVec.size()) _mx = tmpVec.size();
			_chirps.push_back(tmpVec);
		}	
	}	
}

vector<pair<double, double>> Algorithm::add(vector<pair<double, double>> first, vector<pair<double, double>> second)
{
	vector<pair<double, double>> out;
	if (first[0].second < second[0].second)
	{
		out = first;
		for (size_t i = 0; i < second.size(); ++i)
		{
			out.push_back(second[i]);
		}
	}
	else
	{
		out = second;
		for (size_t i = 0; i < first.size(); ++i)
		{
			out.push_back(first[i]);
		}
	}
	sort(out.begin(), out.end(), pairComparison); 
	return out;
}