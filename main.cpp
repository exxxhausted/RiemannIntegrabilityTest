#include <iostream>
#include <cmath>
#include <format>

#include "epsmath.hpp"

using namespace lstu::rit;

struct F {
    double operator () (double x) {
        return 1 / log(x);
    }
};

struct G {
    double operator () (double x) {
        return sin(1/x);
    }
};

struct H {
    double operator () (double x) {
        return fabs(x);
    }
};

int main() {
    double epsilon = 0.01;
    if(auto res = Darboux_kriterium(H(), {-2, 2}, epsilon); !res) {
        std::cout << res.error();
    } else {
        auto& [I, partition] = *res;
        std::cout << std::format("I = ({} +- {})", I, epsilon/2) << std::endl;
        std::cout << partition << std::endl;
    }
    return 0;
}
