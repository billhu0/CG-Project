#include "bezier.h"

#include <vector>
#include <string>

BezierCurve::BezierCurve(int m) { control_points_.resize(m); }

BezierCurve::BezierCurve(std::vector<Vec3f>& control_points, std::vector<float>& weight) {
  control_points_.resize(control_points.size());
  for(int i = 0; i < control_points.size(); ++ i)
  {
    auto control_point = control_points[i];
    auto w = weight[i];
    control_points_[i] = Vec4f(control_point.x() * w, control_point.y() * w, control_point.z() * w, w);
  }
}

BezierCurve::BezierCurve(std::vector<Vec4f>& control_points)
{
  control_points_ = control_points;
}

void BezierCurve::setControlPointAndWeight(int i, Vec3f point, float weight) {
  control_points_[i] = Vec4f(point.x() * weight, point.y() * weight, point.z() * weight, weight);
}

/**
 * TODO: evaluate the point at t with the given control_points
 */
std::pair<Vertex, float> 
BezierCurve::evaluate(std::vector<Vec4f>& control_points, float t) {
  std::vector<Vec4f> P = control_points;
  Vertex res;
  float w;
  for(int j = control_points.size() - 1; j > 0; --j)
  {
    for(int i = 0; i < j; ++ i)
      P[i] = ((float)1.0 - t) * P[i] + t * P[i + 1];
    if(j == 2)
    {
      Vec3f P1(P[1].x() / P[1].w(), P[1].y() / P[1].w(), P[1].z() / P[1].w());
      Vec3f P0(P[0].x() / P[0].w(), P[0].y() / P[0].w(), P[0].z() / P[0].w());
      res.normal = (P1 - P0).normalized();
    }
    if(j == 1)
    {
      w = P[0].w();
      res.position = Vec3f(P[0].x() / w, P[0].y() / w, P[0].z() / w);
    }
  }
  return std::pair(res, w);
}

std::pair<Vertex, float> 
BezierCurve::evaluate(float t) {
  return evaluate(control_points_, t);
}


BezierSurface::BezierSurface(int m, int n, Vec2f range_u, Vec2f range_v) {
  control_points_m_.resize(m);
  for (auto& sub_vec : control_points_m_) {
    sub_vec.resize(n);
  }
  control_points_n_.resize(n);
  for (auto& sub_vec : control_points_n_) {
    sub_vec.resize(m);
  }
  range_u_ = range_u;
  range_v_ = range_v;
}

/**
 * @param[in] i: index (i < m)
 * @param[in] j: index (j < n)
 * @param[in] point: the control point with index i, j
 */
void BezierSurface::setControlPointAndWeight(int i, int j, Vec3f point, float weight) {
  Vec4f wpoint = Vec4f(point.x() * weight, point.y() * weight, point.z() * weight, weight);
  control_points_m_[i][j] = wpoint;
  control_points_n_[j][i] = wpoint;
}

/**
 * TODO: evaluate the point at (u, v) with the given control points(m deriction)
 */
Vertex BezierSurface::evaluate(const std::vector<std::vector<Vec4f>>& control_points, float u, float v) const
{
  std::vector<Vec3f> line;
  std::vector<float> weight;
  Vertex point;
  for(auto sub_vec : control_points)
  {
    BezierCurve subcurve_n(sub_vec);
    auto tmp_controlpoint = subcurve_n.evaluate(v);
    line.push_back(tmp_controlpoint.first.position);
    weight.push_back(tmp_controlpoint.second);
  }
  BezierCurve curve_m(line, weight);
  point = curve_m.evaluate(u).first;
  return point;
}

Vertex BezierSurface::evaluate(float u, float v) const{
  Vertex point_n, point_m, point;
  point_m = evaluate(control_points_m_, u, v);
  point_n = evaluate(control_points_n_, v, u);
  point.position = (point_m.position + point_n.position) * static_cast<float>(0.5);
  point.normal = point_m.normal.cross(point_m.normal);
  return point;
}