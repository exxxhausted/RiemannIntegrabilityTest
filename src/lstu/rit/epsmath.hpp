#ifndef RIT_EPSMATH_HPP
#define RIT_EPSMATH_HPP

#include "Segment.hpp"
#include "Partition.hpp"

#include <expected>
#include <optional>
#include <string>
#include <utility>
#include <map>

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
    double epsilon1 = epsilon/10;
    Partition p(seg, UniformPartitionProperty, 4);
    auto& segs = p.segments();
    std::multimap<double, std::list<Segment>::const_iterator> oscillationMap;
    for(auto it = segs.begin(); it != segs.end(); ++it) {
        auto seg = *it;
        auto Mi = sup_f(f, seg, epsilon1);
        auto mi = inf_f(f, seg, epsilon1);
        if(!Mi || !mi) return std::unexpected("Unintegrable - the gap has been detected.");
        oscillationMap.emplace((((*Mi)-(*mi)) * seg.delta()), it);
    }

    auto S = S_Darboux(f, p, epsilon1);
    auto s = s_Darboux(f, p, epsilon1);
    if(!S) return std::unexpected(S.error());
    if(!s) return std::unexpected(s.error());

    while(true) {
        auto [maxOscillation, segmentIt] = *oscillationMap.rbegin();
        oscillationMap.erase(--oscillationMap.end());
        double oldPoint = segmentIt->a();
        double newPoint = (segmentIt->b() + segmentIt->a()) / 2;
        auto [it, nextIt] = p.addPoint(newPoint);
        auto seg1 = *it;
        auto seg2 = *nextIt;

        auto M1 = sup_f(f, seg1, epsilon1);
        auto m1 = inf_f(f, seg1, epsilon1);
        if(!M1 || !m1) return std::unexpected("Unintegrable - the gap has been detected.");
        oscillationMap.insert( { ((*M1) - (*m1)) * seg1.delta(), it } );

        auto M2 = sup_f(f, seg2, epsilon1);
        auto m2 = inf_f(f, seg2, epsilon1);
        if(!M2 || !m2) return std::unexpected("Unintegrable - the gap has been detected.");
        oscillationMap.insert( { ((*M2) - (*m2)) * seg2.delta(), nextIt } );

        double sum = 0;
        for(const auto& [oscillation, _] : oscillationMap) sum += oscillation;
        if(sum < epsilon) break;
        if(p.min_delta() < epsilon1) return std::unexpected("Unintegrable for this epsilon");
    }

    S = S_Darboux(f, p, epsilon1);
    s = s_Darboux(f, p, epsilon1);
    if(!S) return std::unexpected(S.error());
    if(!s) return std::unexpected(s.error());
    return std::make_pair(((*S) + (*s)) / 2, p);
}

} // lstu::rit

#endif // RIT_EPSMATH_HPP
