#include "bezier.h"

#include <vector>
#include <string>

NURBSPatch::NURBSPatch(int m, int n, Vec2f range_u, Vec2f range_v)
{
  control_points_.resize(m);
  for (auto& sub_vec : control_points_) {
    sub_vec.resize(n);
  }
  degree_u_ = m - 1;
  degree_v_ = n - 1;
  range_u_ = range_u;
  range_v_ = range_v;
  // span_u_ = span_u;
  // span_v_ = span_v;
}

void 
NURBSPatch::setControlPointAndWeight(int i, int j, Vec3f point, float weight)
{
  Vec4f wpoint = Vec4f(point.x() * weight, point.y() * weight, point.z() * weight, weight);
  control_points_[i][j] = wpoint;
}
  
// void 
// NURBSPatch::setParameter(std::vector<float> knots_u, std::vector<float> knots_v)
// {
//   au.resize(degree_u_ + 1);
//   bu.resize(degree_u_ + 1);
//   cu.resize(degree_u_ + 1);
//   du.resize(degree_u_ + 1);
//   for(int i = 1; i <= degree_u_; ++ i)
//   {
//     au[i].resize(degree_u_ + 1);
//     bu[i].resize(degree_u_ + 1);
//     cu[i].resize(degree_u_ + 1);
//     du[i].resize(degree_u_ + 1);
//     for(int j = 0; j <=degree_u_; ++ j)
//     {
//       int k = j + span_u_ - degree_u_;
//       au[i][j] = knots_u[k + i - 1] - knots_u[k];
//       float x = (knots_u[k + i] - knots_u[k + 1]);
//       bu[i][j] = x < EPS ? 0.0f : -1.0f / x;
//       cu[i][j] = -au[i][j] * knots_u[k];
//       du[i][j] = -bu[i][j] * knots_u[k] + i;
//       // printf("%d %d %f %f %f %f\n", i, j, au[i][j], bu[i][j], cu[i][j], du[i][j]);
//     }
//   }
//   av.resize(degree_v_ + 1);
//   bv.resize(degree_v_ + 1);
//   cv.resize(degree_v_ + 1);
//   dv.resize(degree_v_ + 1);
//   for(int i = 1; i <= degree_v_; ++ i)
//   {
//     av[i].resize(degree_v_ + 1);
//     bv[i].resize(degree_v_ + 1);
//     cv[i].resize(degree_v_ + 1);
//     dv[i].resize(degree_v_ + 1);
//     for(int j = 0; j <= degree_v_; ++ j)
//     {
//       int k = j + span_v_ - degree_v_;
//       av[i][j] = knots_v[k + i - 1] - knots_v[k];
//       float x = (knots_v[k + i] - knots_v[k + 1]);
//       bv[i][j] = x < EPS ? 0.0f : -1.0f / x;
//       cv[i][j] = -av[i][j] * knots_v[k];
//       dv[i][j] = -bv[i][j] * knots_v[k] + i;
//     }
//   }
//   return;
// }
  
std::pair<std::vector<float>, std::vector<float>> 
NURBSPatch::EvaluateBasisFunctionDivisionFree_u(float t) const
{
  std::vector<float> N, D;
  N.resize(degree_u_ + 1);
  D.resize(degree_u_ + 1);
  N[0] = 0.0f;
  for(int j = 1; j <= degree_u_; ++ j)
  {
    D[j] = degree_u_ * (au[j][degree_u_] * N[j - 1]);
    N[j] = degree_u_ * (au[j][degree_u_] + cu[j][degree_u_]) * N[j - 1];
    for(int kk = j - 1; kk > 0; -- kk)
    {
      auto k = degree_u_ - j + kk;
      D[kk] = degree_u_ * (au[j][k] * N[kk] - bu[j][k] * N[kk + 1]);
      N[kk] = (au[j][k] * t + cu[j][k]) * N[kk] + (bu[j][k] * t + du[j][k]) * N[kk + 1];
    }
    D[0] = degree_u_ * (-au[j][degree_u_ - j]);
    N[0] = (au[j][degree_u_ - j] * t + cu[j][degree_u_ - j]) * N[j];
  }
  return std::pair(N, D);
}

std::pair<std::vector<float>, std::vector<float>> 
NURBSPatch::EvaluateBasisFunctionDivisionFree_v(float t) const
{
  std::vector<float> N, D;
  N.resize(degree_v_ + 1);
  D.resize(degree_v_ + 1);
  N[0] = 0.0f;
  for(int j = 1; j <= degree_v_; ++ j)
  {
    D[j] = degree_v_ * (av[j][degree_v_] * N[j - 1]);
    N[j] = degree_v_ * (av[j][degree_v_] + cv[j][degree_v_]) * N[j - 1];
    for(int kk = j - 1; kk > 0; -- kk)
    {
      auto k = degree_v_ - j + kk;
      D[kk] = degree_v_ * (av[j][k] * N[kk] - bv[j][k] * N[kk + 1]);
      N[kk] = (av[j][k] * t + cv[j][k]) * N[kk] + (bv[j][k] * t + dv[j][k]) * N[kk + 1];
    }
    D[0] = degree_v_ * (-av[j][degree_v_ - j]);
    N[0] = (av[j][degree_v_ - j] * t + cv[j][degree_v_ - j]) * N[j];
  }
  return std::pair(N, D);
}

std::pair<Vertex, std::pair<Vec3f, Vec3f>> 
NURBSPatch::evaluate(float u, float v) const
{
  // TODO: evaluate the point at (u, v)
  // printf("in %f %f\n", u, v);
  // auto res_u = EvaluateBasisFunctionDivisionFree_u(u);
  // auto res_v = EvaluateBasisFunctionDivisionFree_v(v);
  int span_u_ = find_i(u, knots_u_);
  int span_v_ = find_i(v, knots_v_);
  auto res_u = EvaluateBasisFunctionDirect(u, knots_u_, degree_u_, span_u_);
  auto res_v = EvaluateBasisFunctionDirect(v, knots_v_, degree_v_, span_v_);
  auto Nu = res_u.first;
  auto Du = res_u.second;
  auto Nv = res_v.first;
  auto Dv = res_v.second;

  Vec4f S(0.0f, 0.0f, 0.0f, 0.0f), S_u(0.0f, 0.0f, 0.0f, 0.0f), S_v(0.0f, 0.0f, 0.0f, 0.0f);

  for(int i = 0; i <= degree_u_; ++ i)
    for(int j = 0; j <= degree_v_; ++ j)
    {
      S = S + control_points_[i][j] * Nu[i] * Nv[j];
      S_u = S_u + control_points_[i][j] * Du[i] * Nv[j];
      S_v = S_v + control_points_[i][j] * Nu[i] * Dv[j];
      // printf("%f %f\n", Nu[i], Du[i]);
    }
  Vec3f res_D = Vec3f(S.x(), S.y(), S.z());
  Vec3f res_Du = Vec3f(S_u.x(), S_u.y(), S_u.z());
  Vec3f res_Dv = Vec3f(S_v.x(), S_v.y(), S_v.z());
  Vec3f res_S = res_D / S.w();
  Vec3f res_Su = (res_Du - res_D * S_u.w() / S.w()) / S.w();
  Vec3f res_Sv = (res_Dv - res_D * S_v.w() / S.w()) / S.w();
  Vertex res;
  res.position = res_S;
  res.normal = res_Su.cross(res_Sv).normalized();
  // printf("S: %f %f %f\n", res_S.x(), res_S.y(), res_S.z());
  // printf("Su: %f %f %f\n", res_Su.x(), res_Su.y(), res_Su.z());
  // printf("Sv: %f %f %f\n", res_Sv.x(), res_Sv.y(), res_Sv.z());
  // printf("out\n");
  return std::pair(res, std::pair(res_Su, res_Sv));
}

void
NURBSPatch::setKnots(std::vector<float>& knots_u, std::vector<float>& knots_v)
{
  knots_u_ = knots_u;
  knots_v_ = knots_v;
}

std::pair<std::vector<float>, std::vector<float>> 
NURBSPatch::EvaluateBasisFunctionDirect(float t, std::vector<float> knots, int p, int i) const
{
  //TODO : Compute all non-zero B-spline basis functions and derivatives
  std::vector<float> N, D;
  N.resize(p + 1);
  D.resize(p + 1);
  N[0] = 1.0f;
  std::vector<float> r, l;

  r.resize(p + 1, 0.0f);
  l.resize(p + 1, 0.0f);
  for(int j = 1; j <= p; ++ j)
  {
    float Rn = 0, Rd = 0;
    l[j] = t - knots[i + 1 - j];
    r[j] = knots[i + j] - t;
    for(int k = 0; k <= j - 1; ++ k)
    {
      // printf("%f\n", t);
      
      // printf("a: %f b: %f\n", a, b);
      float Q = N[k] / (r[k + 1] + l[j - k]); 
      // printf("%f\n", Q);
      N[k] = Rn + r[k + 1] * Q;
      Rn = l[j - k] * Q;
      D[k] = p * (Rd - Q);
      Rd = Q;
    }
    N[j] = Rn;
    D[j] = p * Rd;
    // printf("%d %f %f\n", j, N[j], D[j]);
  }
  return std::pair(N, D);
}

int 
NURBSPatch::find_i(float t, std::vector<float> knots) const
{
  for(int i = 0; i < knots.size() - 1; ++ i)
  {
    if(knots[i] <= t && knots[i + 1] > t) return i;
    if(knots[i + 1] == knots[knots.size() - 1] && t == knots[i + 1]) return i;
  }
  return -1;
}