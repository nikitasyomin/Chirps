#pragma once

#include <algorithm>
#include <vector>

using namespace std;

bool pairComparison(pair<double, double> a, pair<double, double> b)
{
	if (a.second != b.second)
		return (a.second < b.second);
	else return (a.first < b.first);
}

bool vectorComparison(vector<pair<double, double>> a, vector<pair<double, double>> b)
{
	if (a[0].second != b[0].second)
		return (a[0].second < b[0].second);
	else return (a[0].first < b[0].first);
}

bool isCorrelated(vector<pair<double, double>> * temp, int first, int last, double R0)
{
	double averageFirst = 0.0;
	double averageSecond = 0.0;
	double R;
	double cov = 0.0;
	double sigmaFirst = 0.0;
	double sigmaSecond = 0.0;

	for (int i = first; i < last; ++i) {
		averageFirst += (*temp)[i].first;
		averageSecond += (*temp)[i].second;
	}

	averageFirst /= (double)(last - first);
	averageSecond /= (double)(last - first);

	for (int i = first; i < last; ++i) {
		cov += ((*temp)[i].first - averageFirst) * ((*temp)[i].second - averageSecond);
		sigmaFirst += ((*temp)[i].first - averageFirst) * ((*temp)[i].first - averageFirst);
		sigmaSecond += ((*temp)[i].second - averageSecond) * ((*temp)[i].second - averageSecond);
	}

	sigmaFirst = sqrt(sigmaFirst);
	sigmaSecond = sqrt(sigmaSecond);
	R = cov / (sigmaFirst * sigmaSecond);

	if (fabs(R) >= R0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

double distance(pair<double, double> x, pair<double, double> y) {
	double d = ((x.first - y.first) * (x.first - y.first) + (x.second - y.second) * (x.second - y.second));
	return sqrt(d);
}

double distance_t(double x, double y) {
	double d = ((x - y) * (x - y) + (x - y) * (x - y));
	return sqrt(d);
}

pair <double, double> linearRegr(vector<pair<double, double>> * temp, int first, int last)
{
	double length = (double) (last - first + 1);
	double angleValue;
	double constant;
	double xSumm = 0.0;
	double ySumm = 0.0;
	double xySumm = 0.0;
	double xxSumm = 0.0;

	for (int i = first; i <= last; ++i)
	{
		ySumm += (*temp)[i].first;
		xSumm += (*temp)[i].second;
		xySumm += (*temp)[i].first * (*temp)[i].second;
		xxSumm += (*temp)[i].second * (*temp)[i].second;
	}

	angleValue = (length * xySumm - xSumm * ySumm) / (length * xxSumm - xSumm * xSumm);
	constant = (ySumm - angleValue * xSumm) / length;
	return {angleValue, constant};
}