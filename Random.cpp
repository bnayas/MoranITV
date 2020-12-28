#include "Random.h"
#include <random>
#include <time.h>
#include <chrono>


std::function<double()> Random::real_unif_dist(double a, double b)
{
    std::default_random_engine generator;
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator.seed(seed);
    auto dis = std::uniform_real_distribution<double>(a, b);
    return (std::bind(dis, generator));
}

std::function<long int()> Random::int_unif_dist(long int  a, long int b)
{
    std::default_random_engine generator;
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator.seed(seed);
    auto dis = std::uniform_int_distribution<>(a, b);
    return (std::bind(dis, generator));
}

double Random::rand()
{
    return(((double) std::rand() / RAND_MAX));
}
int Random::random(int N)
{
    return  ( ((N * (double)std::rand()) - 1) / RAND_MAX);
}

Random::Random()
{
    srand(time(NULL));
}