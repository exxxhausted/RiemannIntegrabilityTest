#ifndef RIT_PARTITION_HPP
#define RIT_PARTITION_HPP

#include <list>
#include <ostream>
#include <set>

#include "Segment.hpp"

namespace lstu::rit {

struct UniformPartitionProperty_t {};
inline static UniformPartitionProperty_t UniformPartitionProperty;

struct EpsilonBasedPartitionProperty_t {};
inline static EpsilonBasedPartitionProperty_t EpsilonBasedPartitionProperty;


class Partition {
public:

    Partition(const Segment& seg);
    Partition(const Segment& seg, UniformPartitionProperty_t, std::size_t N);
    Partition(const Segment& seg, EpsilonBasedPartitionProperty_t, double epsilon);

    const std::list<double>& points() const noexcept;
    const std::list<Segment>& segments() const noexcept;

    auto addPoint(double x) -> std::pair<std::list<Segment>::iterator,
                                         std::list<Segment>::iterator>;
    double lambda() const noexcept;
    double min_delta() const noexcept;

    friend std::ostream& operator <<(std::ostream&, const Partition&);

private:

    std::multiset<double> lengths_;
    std::list<double> points_;
    std::list<Segment> segments_;

};

} // lstu::rit

#endif // RIT_PARTITION_HPP
