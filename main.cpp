#include "curvature.hpp"
#include "cnthd/mesh.hpp"
#include <iostream>

template <typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& arr)
{
    if (arr.empty())
    {
        os << "{}";
    }
    else{
        os << "{";
        for (const auto& a: arr)
        {
            os << a << ", ";
        }
        os << "\b\b}";
    }
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
{
    if (vec.empty())
    {
        os << "{}";
    }
    else{
        os << "{";
        for (const auto& v: vec)
        {
            os << v << ", ";
        }
        os << "\b\b}";
    }
    return os;
}


int main()
{
    // const cnthd::Mesh mesh = cnthd::read_obj("./cnthd/test/data/DDGSpring2016/sphere_small.obj");
    const cnthd::Mesh mesh = cnthd::read_obj("./cnthd/test/data/DDGSpring2016/bunny.obj");
    const auto [K, H] = curvature(mesh);
    std::cout << K << std::endl;
    std::cout << H << std::endl;
    return 0;
}
