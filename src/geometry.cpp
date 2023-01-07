#include "geometry.h"

#include <utility>
#include <iostream>
#include <queue>

constexpr const int debug_flag = 0;

void Mesh::setMaterial(std::shared_ptr<BSDF> &new_bsdf) {
  bsdf = new_bsdf;
}

TriangleMesh::TriangleMesh(std::vector<Vec3f> vertices, std::vector<Vec3f> normals,
                           std::vector<int> v_index, std::vector<int> n_index) :
    vertices(std::move(vertices)),
    normals(std::move(normals)),
    v_indices(std::move(v_index)),
    n_indices(std::move(n_index)) {}

bool TriangleMesh::intersect(Ray &ray, Interaction &interaction) const {
  if (!lbvh.empty()) {
    lbvhHit(interaction, ray);
    // bvhHit(bvh, interaction, ray);
  } else {
    // If you did not implement BVH
    // directly loop through all triangles in the mesh and test intersection for each triangle.
    for (int i = 0; i < v_indices.size() / 3; i++) {
      Vec3i v_idx(v_indices[3 * i], v_indices[3 * i + 1], v_indices[3 * i + 2]);
      Vec3i n_idx(n_indices[3 * i], n_indices[3 * i + 1], n_indices[3 * i + 2]);
      Interaction temp;
      if (intersectOneTriangle(ray, temp, v_idx, n_idx) && (temp.dist < interaction.dist)) {
        interaction = temp;
      }
    }
  }
  return interaction.type != Interaction::Type::NONE;
}

void TriangleMesh::print_triangle_mesh()
{
  for (auto v : vertices) printf("v %f %f %f\n", v.x(), v.y(), v.z());
  for (auto n : normals) printf("n %f %f %f\n", n.x(), n.y(), n.z());
  printf("v_indx : ");
  for (auto index : v_indices) printf ("%d ", index);
  printf("\nn_indx :");
  for (auto index : n_indices) printf ("%d ", index);
}


const AABB TriangleMesh::getTriangle(int pos)
{
  return triangle_AABB[pos];
}

void TriangleMesh::genAABB_for_BVH(BVHNode* now)
{
  now->aabb = getTriangle(now->begin);
  for(int i = 1; i < now->num; ++ i)
  {
    now->aabb = AABB(getTriangle(now->begin + i), now->aabb);
  }
  return;
}

float TriangleMesh::calCost(AABB a, int numa, AABB b, int numb, AABB N)
{
  return a.getSurfaceArea() / N.getSurfaceArea() * numa + b.getSurfaceArea() / N.getSurfaceArea() * numb;
}

int TriangleMesh::getPartitionMethod(BVHNode* now)
{
  int l = now->begin, r = now->begin + now->num - 1;
  int num = now->num;
  int lsize = 1;
  std::vector<AABB> lprefix_AABB, rprefix_AABB;
  lprefix_AABB.push_back(AABB(getTriangle(l)));
  rprefix_AABB.push_back(AABB(getTriangle(r)));
  for(int i = 1; i < num; ++ i)
  {
    lprefix_AABB.push_back(AABB(lprefix_AABB[i - 1], getTriangle(l + i)));
    rprefix_AABB.push_back(AABB(rprefix_AABB[i - 1], getTriangle(r - i)));
  }
  float cost = calCost(lprefix_AABB[0], 1, rprefix_AABB[num - 2], num - 1, now->aabb);
  for(int size = 2; size < num - 1; ++ size)
  {
    float tmpcost = calCost(lprefix_AABB[size - 1], size, rprefix_AABB[num - size - 1], num - size, now->aabb);
    if(cost > tmpcost)
    {
      cost = tmpcost;
      lsize = size;
    }
  }
  return lsize;
}

void TriangleMesh::sort_aabb(int begin, int end, int axis)
{
  if(begin >= end) return;
  // if(debug_flag == 1) printf("%d %d\n", begin, end);
  AABB pivot_aabb = triangle_AABB[begin];
  int pivot_index = triangle_indices[begin];
  int l = begin, r = end, indx = begin;

  for(;r > l;)
  {
    for(;r > l;)
    {
      if(triangle_AABB[r].getCenter()[axis] < pivot_aabb.getCenter()[axis])
      {
        triangle_AABB[l] = triangle_AABB[r];
        triangle_indices[l] = triangle_indices[r];
        indx = r;
        l ++;
        break;
      }
      r --;
    }
    for(; r > l;)
    {
      if(triangle_AABB[l].getCenter()[axis] > pivot_aabb.getCenter()[axis])
      {
        triangle_AABB[r] = triangle_AABB[l];
        triangle_indices[r] = triangle_indices[l];
        indx = l;
        r --;
        break;
      }
      l ++;
    }
  }
  triangle_AABB[indx] = pivot_aabb;
  triangle_indices[indx] = pivot_index;

  sort_aabb(begin, indx - 1, axis);
  sort_aabb(indx + 1, end, axis);
  return;
}

void TriangleMesh::buildBVH_partition(BVHNode* now, int pre_axis)
{
  genAABB_for_BVH(now);
  // printf("%d %d\n", now->begin, now->num);
  // printf("%f %f %f\n", now->aabb.getDist(0), now->aabb.getDist(1), now->aabb.getDist(2));
  if(now->num <= 25) return;
  float dis = now->aabb.getDist(0);
  now->partition_axis = 0;
  for(int i = 1; i < 3; ++ i)
  {
    if(dis < now->aabb.getDist(i))
    {
      dis = now->aabb.getDist(i);
      now->partition_axis = i;
    }
  }

  // sort [begin, begin + num]
  int num = now->num, begin = now->begin;
  // printf("o\n");
  // if(num == 217853 && begin == 0) debug_flag = 1;
  if(pre_axis != now->partition_axis)
      sort_aabb(begin, begin + num - 1, now->partition_axis);
  // printf("f\n");
  int leftnum = getPartitionMethod(now);
  now->left = new BVHNode();
  now->left->num = leftnum;
  now->left->begin = begin;
  now->right = new BVHNode();
  now->right->num = num - leftnum;
  now->right->begin = begin + now->left->num;
  buildBVH_partition(now->left, now->partition_axis);
  buildBVH_partition(now->right, now->partition_axis);
  return;
}

int TriangleMesh::DFS_BVHTree(BVHNode* now)
{
  lbvh.push_back(LBVH());
  int pos = lbvh.size() - 1; 
  if(now->partition_axis == -1)
  {
    lbvh[pos].set(now, -1);
    return 1;
  }
  int lsize = DFS_BVHTree(now->left);
  int rsize = DFS_BVHTree(now->right);
  lbvh[pos].set(now, lsize + 1);
  return lsize + rsize + 1;

}

void TriangleMesh::buildBVH() {
  // TODO: your implementation
  BVHNode *bvh;
  for(int i = 0; i < v_indices.size() / 3; ++ i)
  {
    triangle_AABB.push_back(AABB(vertices[v_indices[3 * i]], vertices[v_indices[3 * i + 1]], vertices[v_indices[3 * i + 2]]));
    triangle_indices.push_back(i);
  }
  bvh = new BVHNode();
  bvh->num = triangle_AABB.size();
  bvh->begin = 0;
  buildBVH_partition(bvh, -1);
  DFS_BVHTree(bvh);
  delete bvh;
  return;
}

bool TriangleMesh::intersectOneTriangle(Ray &ray,
                                        Interaction &interaction,
                                        const Vec3i &v_idx,
                                        const Vec3i &n_idx) const {
  Vec3f v0 = vertices[v_idx[0]];
  Vec3f v1 = vertices[v_idx[1]];
  Vec3f v2 = vertices[v_idx[2]];
  Vec3f v0v1 = v1 - v0;
  Vec3f v0v2 = v2 - v0;
  Vec3f pvec = ray.direction.cross(v0v2);
  float det = v0v1.dot(pvec);

  float invDet = 1.0f / det;

  Vec3f tvec = ray.origin - v0;
  float u = tvec.dot(pvec) * invDet;
  if (u < 0 || u > 1) return false;
  Vec3f qvec = tvec.cross(v0v1);
  float v = ray.direction.dot(qvec) * invDet;
  if (v < 0 || u + v > 1) return false;
  float t = v0v2.dot(qvec) * invDet;
  if (t < ray.t_min || t > ray.t_max) return false;

  interaction.wo = -ray.direction;
  interaction.dist = t;
  interaction.pos = ray(t);
  interaction.normal = (u * normals[n_idx[1]] + v * normals[n_idx[2]]
      + (1 - u - v) * normals[n_idx[0]]).normalized();
  interaction.material = bsdf;
  interaction.type = Interaction::Type::GEOMETRY;
  return true;
}

void TriangleMesh::lbvhHit(Interaction &interaction, Ray &ray) const
{
  std::queue<int> Node_index;
  Node_index.push(0);
  for(; !Node_index.empty();)
  {
    int now_index = Node_index.front();
    Node_index.pop();
    auto now = lbvh[now_index];
    float t_in, t_out;
    if(now.aabb.intersect(ray, &t_in, &t_out))
    {
      if(now.num != -1)
      {
        for(int j = 0; j < now.num; ++ j)
        {
          int i = triangle_indices[now.triangle_offset + j];
          Vec3i v_idx(v_indices[3 * i], v_indices[3 * i + 1], v_indices[3 * i + 2]);
          Vec3i n_idx(n_indices[3 * i], n_indices[3 * i + 1], n_indices[3 * i + 2]);
          Interaction temp;
          if (intersectOneTriangle(ray, temp, v_idx, n_idx) && (temp.dist < interaction.dist)) {
            interaction = temp;
          }
        }
      }
      else
      {
        Node_index.push(now_index + 1);
        Node_index.push(now_index + now.node_offset);
      }
    }
  }
  return;
}

PatchMesh::PatchMesh(std::vector<BezierSurface> &bezier_patches)
{
  patches = bezier_patches;
}

bool
PatchMesh::intersect(Ray &ray, Interaction &interaction) const
{
  return false;
}

void
PatchMesh::buildBVH()
{
  return;
}

void 
PatchMesh::lbvhHit(Interaction &interaction, Ray &ray) const
{
  // TODO: Use the same BVH traversal to traverse the AABB.
  // When reaches a leaf node, call `intersectOnePatch` to determine the nearest intersection.
  std::queue<int> Node_index;
  Node_index.push(0);
  for(; !Node_index.empty();)
  {
    int now_index = Node_index.front();
    Node_index.pop();
    auto now = lbvh[now_index];
    float t_in, t_out;
    if(now.aabb.intersect(ray, &t_in, &t_out))
    {
      if(now.num != -1)
      {
        for(int j = 0; j < now.num; ++ j)
        {
          int i = patch_indices[now.triangle_offset + j];
          Interaction temp;
        
          if (intersectOnePatch(ray, temp, patches[i]) && (temp.dist < interaction.dist)) {
            interaction = temp;
          }
        }
      }
      else
      {
        Node_index.push(now_index + 1);
        Node_index.push(now_index + now.node_offset);
      }
    }
  }
  return;
}

bool
PatchMesh::intersectOnePatch(Ray &ray, Interaction &interaction, const BezierSurface &patch) const
{

  // Write the ray (o+td) as an intersection of two planes.
  // 把 ray (o+td形式) 写成两个平面的交点形式
  Vec3f N1, N2;
  if (std::abs(ray.direction.x()) > std::abs(ray.direction.y()) && std::abs(ray.direction.x()) > std::abs(ray.direction.z())) {
    N1 = Vec3f(ray.direction.y(), -ray.direction.x(), 0);
  } else {
    N1 = Vec3f(0, ray.direction.z(), -ray.direction.y());
  }
  N1.normalize();
  N2 = N1.cross(ray.direction).normalized();
  
  float d1 = -N1.dot(ray.origin);
  float d2 = -N2.dot(ray.origin);

  constexpr float eps = 1e-6;
  constexpr int MAX_ITER = 7;


  // TODO: write Newton Root-Finder here?
  float error_prev = std::numeric_limits<float>::max();
  // u, v are the initial guess of the intersection point
  const float u_initial = 0.5f, v_initial = 0.5f;
  float u = u_initial, v = v_initial;
  for (int iter = 0; iter < MAX_ITER; ++iter) {
    // S = evaluate surface (u, v)
    Vertex S = patch.evaluate(u, v);    
    Vec2f R = Vec2f(N1.dot(S.position) + d1, N2.dot(S.position) + d2);
    float error = std::abs(R.x()) + std::abs(R.y());
    if (error < eps) {

      // Compute the intersection point
      // Update the interaction variables
      interaction.pos = S.position;
      interaction.dist = (S.position - ray.origin).norm();
      interaction.normal = (patch.evaluate(u+eps, v).position - patch.evaluate(u-eps, v).position).cross(patch.evaluate(u, v+eps).position - patch.evaluate(u, v-eps).position).normalized();
      // interaction.material = ? ;
      interaction.wi = ray.direction;
      interaction.wo = -ray.direction;
      interaction.type = Interaction::Type::GEOMETRY;
      return true;
    }
    if (std::abs(error) > error_prev) {
      return false;
    }
    error_prev = std::abs(error);

    // J = compute Jacobian matrix
    // FIXME: 这里是 S_u(u, v), 应该不能直接patch.evaluate?
    Vec2f Fu = {N1.dot(patch.evaluate(u+eps, v).position), N2.dot(patch.evaluate(u+eps, v).position) + d2};
    Vec2f Fv = {N1.dot(patch.evaluate(u, v+eps).position) + d1, N2.dot(patch.evaluate(u, v+eps).position)};
    Mat2f J;
    J.col(0) = Fu;
    J.col(1) = Fv;

    // if J is singular
    if (std::abs(J.determinant()) < eps) {
      u += 0.1 * (u_initial - u) * drand48();
      v += 0.1 * (v_initial - v) * drand48();
    } else {
      u -= J.inverse().col(0).dot(R);
      v -= J.inverse().col(1).dot(R);
    }
  }
  return false;
}


// void TriangleMesh::bvhHit(BVHNode *p, Interaction &interaction,
//                           Ray &ray) const {
//   // TODO: traverse through the bvh and do intersection test efficiently.
//   float t_in, t_out;
//   if(!p->aabb.intersect(ray, &t_in, &t_out)) return;
//   if(p->partition_axis == -1)
//   {
//     for(int j = 0; j < p->num; ++ j)
//     {
//       int i = triangle_indices[p->begin + j];
//       Vec3i v_idx(v_indices[3 * i], v_indices[3 * i + 1], v_indices[3 * i + 2]);
//       Vec3i n_idx(n_indices[3 * i], n_indices[3 * i + 1], n_indices[3 * i + 2]);
//       Interaction temp;
//       if (intersectOneTriangle(ray, temp, v_idx, n_idx) && (temp.dist < interaction.dist)) {
//         interaction = temp;
//       }
//     }
//   }
//   else
//   {
//     bvhHit(p->left, interaction, ray);
//     bvhHit(p->right, interaction, ray);
//   }
//   return;
// }
