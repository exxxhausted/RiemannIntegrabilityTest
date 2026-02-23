#include "Partition.hpp"

#include <algorithm>
#include <format>

namespace lstu::rit {

Partition::Partition(const Segment& seg) {
    points_.push_back(seg.a());
    points_.push_back(seg.b());
    lengths_.insert(seg.b() - seg.a());
    segments_.push_back(seg);
}

Partition::Partition(const Segment& seg, UniformPartitionProperty_t, std::size_t N) {
    double a = seg.a();
    double b = seg.b();

    double step = (b - a) / static_cast<double>(N);

    points_.clear();
    segments_.clear();
    lengths_.clear();

    for (std::size_t i = 0; i <= N; ++i) points_.push_back(a + i * step);

    points_.back() = b;

    auto it = points_.begin();
    auto tmp = it;
    auto next = ++tmp;

    for (; next != points_.end(); ++it, ++next)
    {
        Segment s{*it, *next};
        segments_.push_back(s);
        lengths_.insert(s.delta());
    }
}

Partition::Partition(const Segment& seg, EpsilonBasedPartitionProperty_t, double epsilon) {
    double a = seg.a();
    double b = seg.b();

    points_.clear();
    segments_.clear();
    lengths_.clear();

    for (double x = a; x < b; x += epsilon) points_.push_back(x);

    if (points_.back() != b) points_.push_back(b);

    auto it = points_.begin();
    auto next = std::next(it);

    for (; next != points_.end(); ++it, ++next)
    {
        Segment s{*it, *next};
        segments_.push_back(s);
        lengths_.insert(s.delta());
    }
}

const std::list<double>& Partition::points() const noexcept { return points_; }

const std::list<Segment>& Partition::segments() const noexcept { return segments_; }

auto Partition::addPoint(double x) -> std::pair<std::list<Segment>::iterator,
                                                std::list<Segment>::iterator>
{
    auto pointIt = std::lower_bound(points_.begin(), points_.end(), x);

    if (pointIt == points_.begin() || pointIt == points_.end()) return {};

    auto prevPointIt = std::prev(pointIt);

    auto segIt = segments_.begin();
    for (auto pIt = points_.begin(); std::next(pIt) != pointIt; ++pIt, ++segIt) {}

    Segment oldSeg = *segIt;
    double oldLen = oldSeg.delta();

    auto lenIt = lengths_.find(oldLen);
    if (lenIt != lengths_.end()) lengths_.erase(lenIt);

    segIt = segments_.erase(segIt);

    points_.insert(pointIt, x);

    Segment left  {*prevPointIt, x};
    Segment right {x, *pointIt};

    auto leftIt  = segments_.insert(segIt, left);
    auto rightIt = segments_.insert(std::next(leftIt), right);

    lengths_.insert(left.delta());
    lengths_.insert(right.delta());

    return {leftIt, rightIt};
}

double Partition::lambda() const noexcept { return *lengths_.rbegin(); }

double Partition::min_delta() const noexcept { return *lengths_.begin(); }

std::ostream& operator <<(std::ostream& os, const Partition& p) {
    for(const auto& x_i : p.points_) os << x_i << " ";
    os << std::endl;
    for(const auto& seg_i : p.segments_) os << std::format("[{}, {}]", seg_i.a(), seg_i.b()) << " ";
    return os;
}

} // lstu::rit
