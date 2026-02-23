#include <iostream>
#include <cmath>

#include "epsmath.hpp"

using namespace lstu::rit;

struct F {
    std::optional<double> operator () (double x) {
        if(x == 0) return std::nullopt;
        return 1/(x);
    }
};

struct G {
    std::optional<double> operator () (double x) {
        if(x == 0) return 0;
        return sin(1/x);
    }
};

struct H {
    std::optional<double> operator () (double x) {
        return sqrt(x);
    }
};

struct P {
    std::optional<double> operator () (double x) {
        return x * x;
    }
};

int main()
{
    auto res = Darboux_kriterium(H(), {0, 10}, 1);
    if(!res) std::cout << res.error();
    else {
        auto [I, partition] = *res;
        std::cout << I << std::endl;
        std::cout << partition << std::endl;
    }
    return 0;
}
