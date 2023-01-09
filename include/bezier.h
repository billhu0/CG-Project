#ifndef _BEZIER_H_
#define _BEZIER_H_
#include "core.h"
#include <vector>
#include <utility>

// class BezierCurve {
//  public:
//   std::vector<Vec4f> control_points_;

//   BezierCurve(int m); 
//   BezierCurve(std::vector<Vec3f>& control_points, std::vector<float>& weight); 
//   BezierCurve(std::vector<Vec4f>& control_points);

//   void setControlPointAndWeight(int i, Vec3f point, float weight); 
//   std::pair<Vertex, float>  evaluate(std::vector<Vec4f>& control_points, float t);
//   std::pair<Vertex, float>  evaluate(float t);
// };



// class BezierSurface {
//  public:
//   std::vector<std::vector<Vec4f>> control_points_m_;
//   std::vector<std::vector<Vec4f>> control_points_n_;
//   Vec2f range_u_, range_v_;

//   BezierSurface(int m, int n, Vec2f range_u, Vec2f range_v);
//   void setControlPointAndWeight(int i, int j, Vec3f point, float weight);
//   Vertex evaluate(const std::vector<std::vector<Vec4f>>& control_points, float u, float v) const;
//   std::pair<Vertex, std::pair<Vec3f, Vec3f>> evaluate(float u, float v) const;
// };

class NURBSPatch
{
  public:
    std::vector<std::vector<Vec4f>> control_points_;
    Vec2f range_u_, range_v_;
    std::vector<std::vector<float>> au, av, bu, bv, cu, cv, du, dv;
    int span_u_, span_v_;
    int degree_u_, degree_v_;
    std::vector<float> knots_u_, knots_v_;
    NURBSPatch(int m, int n, int span_u, int span_v, Vec2f range_u, Vec2f range_v);
    void setControlPointAndWeight(int i, int j, Vec3f point, float weight);
    void setParameter(std::vector<float> knots_u, std::vector<float> knots_v);
    void setKnots(std::vector<float>& knots_u, std::vector<float>& knots_v);
    std::pair<std::vector<float>, std::vector<float>> EvaluateBasisFunctionDirect(float t, std::vector<float> knots, int p, int i) const;
    std::pair<std::vector<float>, std::vector<float>> EvaluateBasisFunctionDivisionFree_u(float t) const;
    std::pair<std::vector<float>, std::vector<float>> EvaluateBasisFunctionDivisionFree_v(float t) const;
    std::pair<Vertex, std::pair<Vec3f, Vec3f>> evaluate(float u, float v) const;
    int find_i(float t, std::vector<float> knots) const;
};
#endif