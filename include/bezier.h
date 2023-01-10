#pragma once

#include "core.h"
#include <vector>
#include <utility>

class NURBSPatch {
   public:
    std::vector<std::vector<Vec4f>> control_points_;
    Vec2f range_u_, range_v_;
    std::vector<std::vector<float>> au, av, bu, bv, cu, cv, du, dv;
    int degree_u_, degree_v_;
    std::vector<float> knots_u_, knots_v_;
    NURBSPatch(int m, int n, Vec2f range_u, Vec2f range_v);
    void setControlPointAndWeight(int i, int j, Vec3f point, float weight);
    [[maybe_unused]] void setParameter(std::vector<float> knots_u, std::vector<float> knots_v);
    void setKnots(std::vector<float>& knots_u, std::vector<float>& knots_v);
    [[nodiscard]] static std::pair<std::vector<float>, std::vector<float>> EvaluateBasisFunctionDirect(float t, std::vector<float> knots, int p, int i) ;
    [[maybe_unused]] [[nodiscard]] std::pair<std::vector<float>, std::vector<float>> EvaluateBasisFunctionDivisionFree_u(float t) const;
    [[maybe_unused]] [[nodiscard]] std::pair<std::vector<float>, std::vector<float>> EvaluateBasisFunctionDivisionFree_v(float t) const;
    [[nodiscard]] std::pair<Vertex, std::pair<Vec3f, Vec3f>> evaluate(float u, float v) const;
    [[nodiscard]] static int find_i(float t, std::vector<float> knots) ;
};
