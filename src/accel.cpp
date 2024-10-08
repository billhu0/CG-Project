#include "accel.h"

AABB::AABB(const Vec3f &v1, const Vec3f &v2, const Vec3f &v3) {
    low_bnd = v1.cwiseMin(v2.cwiseMin(v3));
    upper_bnd = v1.cwiseMax(v2.cwiseMax(v3));
}

AABB::AABB(const NURBSPatch &patch) {
    int m = patch.degree_u_;
    int n = patch.degree_v_;
    Vec4f p4 = patch.control_points_[0][0];
    Vec3f p3 = Vec3f(p4.x() / p4.w(), p4.y() / p4.w(), p4.z() / p4.w());
    low_bnd = p3;
    upper_bnd = p3;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            p4 = patch.control_points_[i][j];
            p3 = Vec3f(p4.x() / p4.w(), p4.y() / p4.w(), p4.z() / p4.w());
            low_bnd = low_bnd.cwiseMin(p3);
            upper_bnd = upper_bnd.cwiseMax(p3);
        }
    }
}

AABB::AABB(const AABB &a, const AABB &b) {
    low_bnd = a.low_bnd.cwiseMin(b.low_bnd);
    upper_bnd = a.upper_bnd.cwiseMax(b.upper_bnd);
}

[[maybe_unused]] bool AABB::isOverlap(const AABB &other) const {
    return ((other.low_bnd[0] >= this->low_bnd[0] && other.low_bnd[0] <= this->upper_bnd[0]) ||
            (this->low_bnd[0] >= other.low_bnd[0] && this->low_bnd[0] <= other.upper_bnd[0])) &&
           ((other.low_bnd[1] >= this->low_bnd[1] && other.low_bnd[1] <= this->upper_bnd[1]) ||
            (this->low_bnd[1] >= other.low_bnd[1] && this->low_bnd[1] <= other.upper_bnd[1])) &&
           ((other.low_bnd[2] >= this->low_bnd[2] && other.low_bnd[2] <= this->upper_bnd[2]) ||
            (this->low_bnd[2] >= other.low_bnd[2] && this->low_bnd[2] <= other.upper_bnd[2]));
}

bool AABB::intersect(const Ray &ray, float *t_in, float *t_out) {
    // intersection test for bounding box
    // ray distance for two intersection points are returned by pointers.
    float dir_frac_x = (ray.direction[0] == 0.0) ? 1.0e32f : 1.0f / ray.direction[0];
    float dir_frac_y = (ray.direction[1] == 0.0) ? 1.0e32f : 1.0f / ray.direction[1];
    float dir_frac_z = (ray.direction[2] == 0.0) ? 1.0e32f : 1.0f / ray.direction[2];

    float tx1 = (low_bnd[0] - ray.origin[0]) * dir_frac_x;
    float tx2 = (upper_bnd[0] - ray.origin[0]) * dir_frac_x;
    float ty1 = (low_bnd[1] - ray.origin[1]) * dir_frac_y;
    float ty2 = (upper_bnd[1] - ray.origin[1]) * dir_frac_y;
    float tz1 = (low_bnd[2] - ray.origin[2]) * dir_frac_z;
    float tz2 = (upper_bnd[2] - ray.origin[2]) * dir_frac_z;

    *t_in = std::max(std::max(std::min(tx1, tx2), std::min(ty1, ty2)), std::min(tz1, tz2));
    *t_out = std::min(std::min(std::max(tx1, tx2), std::max(ty1, ty2)), std::max(tz1, tz2));

    // When t_out < 0 and the ray is intersecting with AABB, the whole AABB is behind us
    *t_in = std::max(*t_in, ray.t_min);
    *t_out = std::min(*t_out, ray.t_max);
    if (*t_out < 0) return false;

    return *t_out >= *t_in;
}

float AABB::getSurfaceArea() const {
    float x = upper_bnd.x() - low_bnd.x(), y = upper_bnd.y() - low_bnd.y(), z = upper_bnd.z() - low_bnd.z();
    return 2.0f * (x * y + x * z + y * z);
}
