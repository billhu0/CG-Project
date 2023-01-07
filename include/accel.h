#ifndef ACCEL_H_
#define ACCEL_H_

#include "core.h"
#include "ray.h"
#include "bezier.h"

struct AABB {
  // the minimum and maximum coordinate for the AABB
  Vec3f low_bnd;
  Vec3f upper_bnd;
  /// test intersection with given ray.
  /// ray distance of entrance and exit point are recorded in t_in and t_out
  AABB() : low_bnd(0, 0, 0), upper_bnd(0, 0, 0) {}
  AABB(Vec3f low, Vec3f upper) : low_bnd(low), upper_bnd(upper) {}
  /// construct an AABB from three vertices of a triangle.
  AABB(const Vec3f &v1, const Vec3f &v2, const Vec3f &v3);
  /// construct an AABB from the bezier surface
  AABB(const BezierSurface &patch); 
  /// Construct AABB by merging two AABBs
  AABB(const AABB &a, const AABB &b);
  bool intersect(const Ray &ray, float *t_in, float *t_out);
  /// Get the AABB center
  [[nodiscard]] Vec3f getCenter() const { return (low_bnd + upper_bnd) / 2; }
  /// Get the length of a specified side on the AABB
  [[nodiscard]] float getDist(int dim) const { return upper_bnd[dim] - low_bnd[dim]; }
  /// Check whether the AABB is overlapping with another AABB
  [[nodiscard]] bool isOverlap(const AABB &other) const;
  [[nodiscard]] float getSurfaceArea() const;
};

struct BVHNode {
  BVHNode *left;
  BVHNode *right;
  // bounding box of current node.
  AABB aabb;
  int partition_axis; 
  // index of triangles in current BVH leaf node.
  // std::vector<int> triangles;
  int num, begin;

  BVHNode() : aabb(), left(nullptr), right(nullptr), partition_axis(-1), num(0), begin(0) {}
  ~BVHNode()
  {
    if(left) delete left;
    if(right) delete right;
  }
};

struct LBVH{
  AABB aabb;
  int num; //positiv : leaf_size, -1 : not leaf node
  union {
    int node_offset;
    int triangle_offset;
  } ;


  void set(BVHNode *node, int in_offset)
  {
    if(in_offset != -1)
    {
      num = -1;
      aabb = node->aabb;
      node_offset = in_offset;
    }
    else
    {
      num = node->num;
      aabb = node->aabb;
      triangle_offset = node->begin;
    }
    return;
  }
};
// You may need to add your code for BVH construction here.

#endif //ACCEL_H_