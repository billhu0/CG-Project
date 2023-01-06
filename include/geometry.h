#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include "core.h"
#include "ray.h"
#include "interaction.h"
#include "bsdf.h"
#include "accel.h"

#include <vector>
#include <optional>

class TriangleMesh {
 public:
  TriangleMesh() = default;
  TriangleMesh(std::vector<Vec3f> vertices,
               std::vector<Vec3f> normals,
               std::vector<int> v_index,
               std::vector<int> n_index);
  ~TriangleMesh(){}
  bool intersect(Ray &ray, Interaction &interaction) const;
  void setMaterial(std::shared_ptr<BSDF> &new_bsdf);
  const AABB getTriangle(int pos);
  void buildBVH();
  void print_triangle_mesh();
 private:
  void buildBVH_partition(BVHNode* now, int pre_axis);
  int getPartitionMethod(BVHNode* now);
  float calCost(AABB a, int numa, AABB b, int numb, AABB N);
  int DFS_BVHTree(BVHNode* now);
  void genAABB_for_BVH(BVHNode* now);
  void sort_aabb(int begin, int end, int axis);
  bool intersectOneTriangle(Ray &ray, Interaction &interaction, const Vec3i& v_idx, const Vec3i& n_idx) const;
  void lbvhHit(Interaction &Interaction, Ray &ray) const;
  // void bvhHit(BVHNode *p, Interaction &interaction, Ray &ray) const;
  
  std::shared_ptr<BSDF> bsdf;
  // BVHNode* bvh;
  

  std::vector<Vec3f> vertices;
  std::vector<Vec3f> normals;
  std::vector<int> v_indices;
  std::vector<int> n_indices;
  std::vector<AABB> triangle_AABB;
  std::vector<int> triangle_indices;
  std::vector<LBVH> lbvh;
};

#endif // GEOMETRY_H_
