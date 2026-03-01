#ifndef RIT_EPSMATH_HPP
#define RIT_EPSMATH_HPP

#include "Segment.hpp"
#include "Partition.hpp"

#include <expected>
#include <optional>
#include <string>
#include <utility>
#include <map>
#include <cmath>
#include <limits>

namespace lstu::rit {

template<typename F>
double inf_f(F f, const Segment& seg, double epsilon) {
    Partition p(seg, EpsilonBasedPartitionProperty, epsilon);

    double min = std::numeric_limits<double>::quiet_NaN();

    for (const auto& x : p.points()) {
        double value = f(x);
        if (!std::isnan(value)) {
            min = value;
            break;
        }
    }

    if (std::isnan(min)) return std::numeric_limits<double>::quiet_NaN();

    for (const auto& x : p.points()) {
        double value = f(x);
        if (!std::isnan(value) && value < min) min = value;
    }

    return min;
}

template<typename F>
double sup_f(F f, const Segment& seg, double epsilon) {
    Partition p(seg, EpsilonBasedPartitionProperty, epsilon);

    double max = std::numeric_limits<double>::quiet_NaN();

    for (const auto& x : p.points()) {
        double value = f(x);
        if (!std::isnan(value)) {
            max = value;
            break;
        }
    }

    if (std::isnan(max)) return std::numeric_limits<double>::quiet_NaN();

    for (const auto& x : p.points()) {
        double value = f(x);
        if (!std::isnan(value) && value > max) max = value;
    }

    return max;
}

template<typename F>
std::expected<double, std::string> S_Darboux(F f, const Partition& p, double epsilon) {
    double sum = 0;
    for(const auto& seg : p.segments()) {
        double Mi = sup_f(f, seg, epsilon);
        if(std::isnan(Mi)) return std::unexpected("Darboux sum doesent exists - the gap has been detected.");
        if(std::isinf(Mi)) return std::unexpected("Darboux sum doesent exists - the function is not limited.");
        sum += Mi * seg.delta();
    }
    return sum;
}

template<typename F>
std::expected<double, std::string> s_Darboux(F f, const Partition& p, double epsilon) {
    double sum = 0;
    for(const auto& seg : p.segments()) {
        double mi = inf_f(f, seg, epsilon);
        if(std::isnan(mi)) return std::unexpected("Darboux sum doesent exists - the gap has been detected.");
        if(std::isinf(mi)) return std::unexpected("Darboux sum doesent exists - the function is not limited.");
        sum += mi * seg.delta();
    }
    return sum;
}

template<typename F>
auto Darboux_kriterium(F f, const Segment& seg, double epsilon) ->
    std::expected<std::pair<double, Partition>, std::string>
{
    Partition p(seg, UniformPartitionProperty, 2);
    auto& segs = p.segments();
    std::multimap<double, std::list<Segment>::const_iterator> oscillationMap;

    double oscillationSum = 0;
    for(auto it = segs.begin(); it != segs.end(); ++it) {
        auto seg = *it;
        double Mi = sup_f(f, seg, seg.delta() / 10);
        if(std::isnan(Mi)) return std::unexpected("Not Riemann integrable - the gap has been detected.");
        if(std::isinf(Mi)) return std::unexpected("Not Riemann integrable - the function is not limited.");
        double mi = inf_f(f, seg, seg.delta() / 10);
        if(std::isnan(mi)) return std::unexpected("Not Riemann integrable - the gap has been detected.");
        if(std::isinf(mi)) return std::unexpected("Not Riemann integrable - the function is not limited.");
        double value = (Mi - mi) * seg.delta();
        oscillationMap.emplace(value, it);
        oscillationSum += value;
    }

    while(true) {
        if (oscillationSum < epsilon) break;
        if (p.min_delta() < std::numeric_limits<double>::epsilon())
            return std::unexpected("Not Riemann integrable - machine precision limit");

        auto [maxOscillation, segmentIt] = *oscillationMap.rbegin();
        oscillationMap.erase(--oscillationMap.end());
        double oldPoint = segmentIt->a();
        double newPoint = (segmentIt->b() + segmentIt->a()) / 2;
        auto [it, nextIt] = p.addPoint(newPoint);
        oscillationSum -= maxOscillation;
        auto seg1 = *it;
        auto seg2 = *nextIt;

        double M1 = sup_f(f, seg1, seg1.delta() / 10);
        if(std::isnan(M1)) return std::unexpected("Not Riemann integrable - the gap has been detected.");
        if(std::isinf(M1)) return std::unexpected("Not Riemann integrable - the function is not limited.");
        double m1 = inf_f(f, seg1, seg1.delta() / 10);
        if(std::isnan(m1)) return std::unexpected("Not Riemann integrable - the gap has been detected.");
        if(std::isinf(m1)) return std::unexpected("Not Riemann integrable - the function is not limited.");
        double value1 = (M1 - m1) * seg1.delta();
        oscillationMap.insert( { value1, it } );
        oscillationSum += value1;

        double M2 = sup_f(f, seg2, seg2.delta() / 10);
        if(std::isnan(M2)) return std::unexpected("Not Riemann integrable - the gap has been detected.");
        if(std::isinf(M2)) return std::unexpected("Not Riemann integrable - the function is not limited.");
        double m2 = inf_f(f, seg2, seg2.delta() / 10);
        if(std::isnan(m2)) return std::unexpected("Not Riemann integrable - the gap has been detected.");
        if(std::isinf(m2)) return std::unexpected("Not Riemann integrable - the function is not limited.");
        double value2 = (M2 - m2) * seg2.delta();
        oscillationMap.insert( { value2, nextIt } );
        oscillationSum += value2;
    }

    auto S = S_Darboux(f, p, epsilon);
    auto s = s_Darboux(f, p, epsilon);
    if(!S) return std::unexpected(S.error());
    if(!s) return std::unexpected(s.error());
    return std::make_pair(((*S) + (*s)) / 2, p);
}

} // lstu::rit

#endif // RIT_EPSMATH_HPP
