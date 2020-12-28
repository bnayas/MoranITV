#pragma once
#include <functional>

class Random
{
public:
	Random();
	std::function<double()> real_unif_dist(double a, double b);
	std::function<long int()> int_unif_dist(long int a, long int b);
	double rand();
	int random(int N);
};

