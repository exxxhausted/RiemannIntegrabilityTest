#include <iostream>
#include <cmath>

#include "epsmath.hpp"

using namespace lstu::rit;

struct F {
    std::optional<double> operator () (double x) {
        if(x == 0) return 0;
        else return sin(1/x);
    }
};

int main()
{
    auto res = Darboux_kriterium(F(), {-3.14, 3.14}, 0.001);
    if(!res) std::cout << res.error();
    else {
        auto [I, partition] = *res;
        std::cout << I << std::endl;
        std::cout << partition << std::endl;
    }
    return 0;
}
