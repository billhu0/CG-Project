#ifndef CORE_H_
#define CORE_H_

#include <Eigen/Core>
#include <Eigen/Dense>
#include <memory>

#define UNIMPLEMENTED

template <typename T, int size>
using Vec = Eigen::Matrix<T, size, 1>;

using Vec3f = Eigen::Vector3f;
using Vec3i = Eigen::Vector3i;
using Vec2f = Eigen::Vector2f;
using Vec2i = Eigen::Vector2i;
using Vec4f = Eigen::Vector4f;
using Mat3f = Eigen::Matrix3f;
using Mat2f = Eigen::Matrix2f;

constexpr float RAY_DEFAULT_MIN = 1e-5;
constexpr float RAY_DEFAULT_MAX = 1e7;
constexpr float PI = 3.141592653579f;
constexpr float INV_PI = 0.31830988618379067154;
constexpr float EPS = 1e-5;
constexpr bool ENABLE_NURBS = true;
class BSDF;
class Sampler;

struct Vertex {
    Vec3f position;
    Vec3f normal;
};

#endif  // CORE_H_
