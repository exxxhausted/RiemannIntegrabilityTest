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

namespace lstu::rit {

template<typename F>
std::optional<double> inf_f(F f, const Segment& seg, double epsilon) {
    Partition p(seg, EpsilonBasedPartitionProperty, epsilon);

    std::optional<double> min = std::nullopt;

    for (const auto& x : p.points()) {
        auto value = f(x);
        if (value) {
            min = value;
            break;
        }
    }

    if (!min)
        return std::nullopt;

    for (const auto& x : p.points()) {
        auto value = f(x);
        if (value && *value < *min)
            min = value;
    }

    return min;
}

template<typename F>
std::optional<double> sup_f(F f, const Segment& seg, double epsilon) {
    Partition p(seg, EpsilonBasedPartitionProperty, epsilon);

    std::optional<double> max = std::nullopt;

    for (const auto& x : p.points()) {
        auto value = f(x);
        if (value) {
            max = value;
            break;
        }
    }

    if (!max)
        return std::nullopt;

    for (const auto& x : p.points()) {
        auto value = f(x);
        if (value && *value > *max)
            max = value;
    }

    return max;
}

template<typename F>
std::expected<double, std::string> S_Darboux(F f, const Partition& p, double epsilon) {
    double sum = 0;
    for(const auto& seg : p.segments()) {
        auto Mi = sup_f(f, seg, epsilon);
        if(!Mi) return std::unexpected("Unintegrable - the gap has been detected.");
        sum += *Mi * seg.delta();
    }
    return sum;
}

template<typename F>
std::expected<double, std::string> s_Darboux(F f, const Partition& p, double epsilon) {
    double sum = 0;
    for(const auto& seg : p.segments()) {
        auto mi = inf_f(f, seg, epsilon);
        if(!mi) return std::unexpected("Unintegrable - the gap has been detected.");
        sum += *mi * seg.delta();
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
    double prevOscilationSum = 0;
    for(auto it = segs.begin(); it != segs.end(); ++it) {
        auto seg = *it;
        auto Mi = sup_f(f, seg, epsilon);
        auto mi = inf_f(f, seg, epsilon);
        if(!Mi || !mi) return std::unexpected("Unintegrable - the gap has been detected.");
        double value = ((*Mi)-(*mi)) * seg.delta();
        oscillationMap.emplace(value, it);
        oscillationSum += value;
    }

    auto S = S_Darboux(f, p, epsilon);
    auto s = s_Darboux(f, p, epsilon);
    if(!S) return std::unexpected(S.error());
    if(!s) return std::unexpected(s.error());

    while(true) {

        if(oscillationSum < epsilon) break;
        if (p.min_delta() < std::numeric_limits<double>::epsilon())
            return std::unexpected("Unintegrable - machine precision limit");

        prevOscilationSum = oscillationSum;
        auto [maxOscillation, segmentIt] = *oscillationMap.rbegin();
        oscillationMap.erase(--oscillationMap.end());
        double oldPoint = segmentIt->a();
        double newPoint = (segmentIt->b() + segmentIt->a()) / 2;
        auto [it, nextIt] = p.addPoint(newPoint);
        oscillationSum -= maxOscillation;
        auto seg1 = *it;
        auto seg2 = *nextIt;

        auto M1 = sup_f(f, seg1, seg1.delta() / 100);
        auto m1 = inf_f(f, seg1, seg1.delta() / 100);
        if(!M1 || !m1) return std::unexpected("Unintegrable - the gap has been detected.");
        double value1 = ((*M1) - (*m1)) * seg1.delta();
        oscillationMap.insert( { value1, it } );
        oscillationSum += value1;

        auto M2 = sup_f(f, seg2, seg2.delta() / 100);
        auto m2 = inf_f(f, seg2, seg2.delta() / 100);
        if(!M2 || !m2) return std::unexpected("Unintegrable - the gap has been detected.");
        double value2 = ((*M2) - (*m2)) * seg2.delta();
        oscillationMap.insert( { value2, nextIt } );
        oscillationSum += value2;
    }

    S = S_Darboux(f, p, epsilon);
    s = s_Darboux(f, p, epsilon);
    if(!S) return std::unexpected(S.error());
    if(!s) return std::unexpected(s.error());
    return std::make_pair(((*S) + (*s)) / 2, p);
}

} // lstu::rit

#endif // RIT_EPSMATH_HPP
