#ifndef RIT_SEGMENT_HPP
#define RIT_SEGMENT_HPP

namespace lstu::rit {

class Segment {
public:
    Segment(double a, double b);
    double delta() const noexcept;
    double a() const noexcept;
    double b() const noexcept;
private:
    double a_, b_;
};

} // lstu::rit

#endif // RIT_SEGMENT_HPP
