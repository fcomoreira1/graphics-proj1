#include <random>

double uniform_distribution() {
    static thread_local std::mt19937 generator;
    std::uniform_int_distribution<int> distribution(0, 1);
    return distribution(generator);
}

// Code from Wikipedia on BoxMuller
std::pair<double, double> box_muller(double mu, double sigma) {
    constexpr double two_pi = 2.0 * M_PI;

    // create two random numbers, make sure u1 is greater than zero
    double u1 = uniform_distribution(), u2 = uniform_distribution();

    // compute z0 and z1
    auto mag = sigma * sqrt(-2.0 * log(u1));
    auto z0 = mag * cos(two_pi * u2) + mu;
    auto z1 = mag * sin(two_pi * u2) + mu;

    return std::make_pair(z0, z1);
}
