#include "Segment.hpp"

#include <stdexcept>

namespace lstu::rit {

Segment::Segment(double a, double b) {
    if(a > b) throw std::invalid_argument("b must be greater than a!");
    a_ = a;
    b_ = b;
}

double Segment::delta() const noexcept { return b_ - a_; }

double Segment::a() const noexcept { return a_; }

double Segment::b() const noexcept { return b_; }

} // lstu::rit
