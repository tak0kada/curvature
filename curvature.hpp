#pragma once

#include <array>
#include <vector>
// #include <iostream>
#include <cstddef>
#include <cmath>
#include <boost/math/constants/constants.hpp>
#ifdef _OPENMP
    #include <omp.h>
#endif
#include "cnthd/mesh.hpp"

constexpr double PI = boost::math::constants::pi<double>();

// MEYER, Mark, et al. Discrete differential-geometry operators for triangulated 2-manifolds. In: Visualization and mathematics III. Springer, Berlin, Heidelberg, 2003. p. 35-57.
// std::pair<GaussianCurvature, MeanCurvature>
std::pair<std::vector<double>, std::vector<double>> curvature(const cnthd::Mesh& mesh)
{
    std::vector<double> K(mesh.nV);
    std::vector<double> H(mesh.nV);

    std::vector<double> face_area(mesh.nF);
    for (std::size_t fi = 0; fi < mesh.nF; ++fi)
    {
        face_area[fi] = mesh.face[fi].area();
    }

#pragma omp parallel for
    for (std::size_t vi = 0; vi < mesh.nV; ++vi)
    {
        std::vector<cnthd::HalfEdge*> hes;
        cnthd::HalfEdge* tmp = mesh.vertex[vi].he_out;
        while(true)
        {
            hes.push_back(tmp);
            tmp = tmp->prev->pair;
            if (tmp->idx == mesh.vertex[vi].he_out->idx)
            {
                break;
            }
        }

        // calculate cotangent ------------------------------------------------
        std::vector<double> cota(hes.size());
        std::vector<double> cotb(hes.size());
        for (int hi = 0; hi < (int)hes.size(); ++hi)
        {
            // cos * |u| * |v|
            // sin * |u| * |v|
            double ca = (*hes[hi]->from - *hes[hi]->to) * (*hes[hi]->next->to - *hes[hi]->next->from);
            double sa = ((*hes[hi]->from - *hes[hi]->to) ^ (*hes[hi]->next->to - *hes[hi]->next->from)).length();
            cota[(hi + 1) % hes.size()] =  ca / sa;

            double cb = (*hes[hi]->from - *hes[hi]->to) * (*hes[hi]->pair->prev->from - *hes[hi]->pair->prev->to);
            double sb = ((*hes[hi]->from - *hes[hi]->to) ^ (*hes[hi]->pair->prev->from - *hes[hi]->pair->prev->to)).length();
            cotb[(hi - 1 + hes.size()) % hes.size()] = cb / sb;
        }


        // calculate area -----------------------------------------------------
        double Amixed{0};
        for (int hi = 0; hi < (int)hes.size(); ++hi)
        {
            // if non-obtuse
            if ((*hes[hi]->pair->next->to - *hes[hi]->pair->next->from) * (*hes[hi]->to - *hes[hi]->from) >= 0
                && ((*hes[hi]->to - *hes[hi]->from) * (*hes[hi]->prev->from - *hes[hi]->prev->to)) >= 0)
            {
                Amixed += 0.125 * hes[hi]->edge->norm() * cota[(hi + 1) % hes.size()]
                       + 0.125 * hes[(hi + 1) % hes.size()]->edge->norm() * cotb[hi];
                continue;
            }

            // obtuse
            // if angle(vi) > pi/2
            if (cota[(hi + 1) % hes.size()] > 0 && cotb[hi] > 0)
            {
                Amixed += 0.5 * face_area[hes[hi]->face->idx];
            }
            else
            {
                Amixed += 0.25 * face_area[hes[hi]->face->idx];
            }
        }

        // Mean Curvature
        cnthd::Vector h{0, 0, 0};
        // face normal
        cnthd::Vector n{0, 0, 0};
        for (int hi = 0; hi < (int)hes.size(); ++hi)
        {
            h += (*hes[hi]->to - *hes[hi]->from) * (cota[hi] + cotb[hi]);
            n -= (*hes[hi]->to - *hes[hi]->from) ^ (*hes[hi]->prev->from - *hes[hi]->prev->to);
        }
        int sgn = n*h > 0 ? 1 : -1;
        H[vi] = sgn * h.length() * 0.25 / Amixed;

        // Gaussian Curvature
        double angle_defect{2 * PI};
        for (int hi = 0; hi < (int)hes.size(); ++hi)
        {
            angle_defect -= PI - std::atan2(1, cota[hi]) - std::atan2(1, cotb[hi]);
        }
        K[vi] = angle_defect / Amixed;
    }

    return {K, H};
}
