#include "core.h"
#include "geometry.h"

#include <utility>
#include <iostream>
#include <queue>

void Mesh::setMaterial(std::shared_ptr<BSDF> &new_bsdf) {
    bsdf = new_bsdf;
}

TriangleMesh::TriangleMesh(std::vector<Vec3f> vertices, std::vector<Vec3f> normals, std::vector<int> v_index,
                           std::vector<int> n_index)
    : vertices(std::move(vertices)),
      normals(std::move(normals)),
      v_indices(std::move(v_index)),
      n_indices(std::move(n_index)) {
}

bool TriangleMesh::intersect(Ray &ray, Interaction &interaction) const {
    if (!lbvh.empty()) {
        lbvhHit(interaction, ray);
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

void TriangleMesh::print_triangle_mesh() {
    for (auto v : vertices) printf("v %f %f %f\n", v.x(), v.y(), v.z());
    for (auto n : normals) printf("n %f %f %f\n", n.x(), n.y(), n.z());
    printf("v_indx : ");
    for (auto index : v_indices) printf("%d ", index);
    printf("\nn_indx :");
    for (auto index : n_indices) printf("%d ", index);
}

AABB TriangleMesh::getTriangle(int pos) {
    return triangle_AABB[pos];
}

void TriangleMesh::genAABB_for_BVH(BVHNode *now) {
    now->aabb = getTriangle(now->begin);
    for (int i = 1; i < now->num; ++i) {
        now->aabb = AABB(getTriangle(now->begin + i), now->aabb);
    }
}

float TriangleMesh::calCost(const AABB& a, int numa, const AABB& b, int numb, const AABB& N) {
    return a.getSurfaceArea() / N.getSurfaceArea() * numa + b.getSurfaceArea() / N.getSurfaceArea() * numb;
}

int TriangleMesh::getPartitionMethod(BVHNode *now) {
    int l = now->begin, r = now->begin + now->num - 1;
    int num = now->num;
    int lsize = 1;
    std::vector<AABB> lprefix_AABB, rprefix_AABB;
    lprefix_AABB.push_back(AABB(getTriangle(l)));
    rprefix_AABB.push_back(AABB(getTriangle(r)));
    for (int i = 1; i < num; ++i) {
        lprefix_AABB.emplace_back(lprefix_AABB[i - 1], getTriangle(l + i));
        rprefix_AABB.emplace_back(rprefix_AABB[i - 1], getTriangle(r - i));
    }
    float cost = calCost(lprefix_AABB[0], 1, rprefix_AABB[num - 2], num - 1, now->aabb);
    for (int size = 2; size < num - 1; ++size) {
        float tmpcost = calCost(lprefix_AABB[size - 1], size, rprefix_AABB[num - size - 1], num - size, now->aabb);
        if (cost > tmpcost) {
            cost = tmpcost;
            lsize = size;
        }
    }
    return lsize;
}

void TriangleMesh::sort_aabb(int begin, int end, int axis) {
    if (begin >= end) return;
    // if(debug_flag == 1) printf("%d %d\n", begin, end);
    AABB pivot_aabb = triangle_AABB[begin];
    int pivot_index = triangle_indices[begin];
    int l = begin, r = end, indx = begin;

    for (; r > l;) {
        for (; r > l;) {
            if (triangle_AABB[r].getCenter()[axis] < pivot_aabb.getCenter()[axis]) {
                triangle_AABB[l] = triangle_AABB[r];
                triangle_indices[l] = triangle_indices[r];
                indx = r;
                l++;
                break;
            }
            r--;
        }
        for (; r > l;) {
            if (triangle_AABB[l].getCenter()[axis] > pivot_aabb.getCenter()[axis]) {
                triangle_AABB[r] = triangle_AABB[l];
                triangle_indices[r] = triangle_indices[l];
                indx = l;
                r--;
                break;
            }
            l++;
        }
    }
    triangle_AABB[indx] = pivot_aabb;
    triangle_indices[indx] = pivot_index;

    sort_aabb(begin, indx - 1, axis);
    sort_aabb(indx + 1, end, axis);
}

void TriangleMesh::buildBVH_partition(BVHNode *now, int pre_axis) {
    genAABB_for_BVH(now);
    // printf("%d %d\n", now->begin, now->num);
    // printf("%f %f %f\n", now->aabb.getDist(0), now->aabb.getDist(1), now->aabb.getDist(2));
    if (now->num <= 25) return;
    float dis = now->aabb.getDist(0);
    now->partition_axis = 0;
    for (int i = 1; i < 3; ++i) {
        if (dis < now->aabb.getDist(i)) {
            dis = now->aabb.getDist(i);
            now->partition_axis = i;
        }
    }

    // sort [begin, begin + num]
    int num = now->num, begin = now->begin;
    // printf("o\n");
    // if(num == 217853 && begin == 0) debug_flag = 1;
    if (pre_axis != now->partition_axis) sort_aabb(begin, begin + num - 1, now->partition_axis);
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
}

int TriangleMesh::DFS_BVHTree(BVHNode *now) {
    lbvh.emplace_back();
    int pos = lbvh.size() - 1;
    if (now->partition_axis == -1) {
        lbvh[pos].set(now, -1);
        return 1;
    }
    int lsize = DFS_BVHTree(now->left);
    int rsize = DFS_BVHTree(now->right);
    lbvh[pos].set(now, lsize + 1);
    return lsize + rsize + 1;
}

void TriangleMesh::buildBVH() {
    BVHNode *bvh;
    for (int i = 0; i < v_indices.size() / 3; ++i) {
        triangle_AABB.emplace_back(vertices[v_indices[3 * i]], vertices[v_indices[3 * i + 1]], vertices[v_indices[3 * i + 2]]);
        triangle_indices.push_back(i);
    }
    bvh = new BVHNode();
    bvh->num = (int) triangle_AABB.size();
    bvh->begin = 0;
    buildBVH_partition(bvh, -1);
    DFS_BVHTree(bvh);
    delete bvh;
}

bool TriangleMesh::intersectOneTriangle(Ray &ray, Interaction &interaction, const Vec3i &v_idx, const Vec3i &n_idx) const {
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
    interaction.normal = (u * normals[n_idx[1]] + v * normals[n_idx[2]] + (1 - u - v) * normals[n_idx[0]]).normalized();
    interaction.material = bsdf;
    interaction.type = Interaction::Type::GEOMETRY;
    return true;
}

void TriangleMesh::lbvhHit(Interaction &interaction, Ray &ray) const {
    std::queue<int> Node_index;
    Node_index.push(0);
    for (; !Node_index.empty();) {
        int now_index = Node_index.front();
        Node_index.pop();
        auto now = lbvh[now_index];
        float t_in, t_out;
        if (now.aabb.intersect(ray, &t_in, &t_out)) {
            if (now.num != -1) {
                for (int j = 0; j < now.num; ++j) {
                    int i = triangle_indices[now.object_offset + j];
                    Vec3i v_idx(v_indices[3 * i], v_indices[3 * i + 1], v_indices[3 * i + 2]);
                    Vec3i n_idx(n_indices[3 * i], n_indices[3 * i + 1], n_indices[3 * i + 2]);
                    Interaction temp;
                    if (intersectOneTriangle(ray, temp, v_idx, n_idx) && (temp.dist < interaction.dist)) {
                        interaction = temp;
                    }
                }
            } else {
                Node_index.push(now_index + 1);
                Node_index.push(now_index + now.node_offset);
            }
        }
    }
}

PatchMesh::PatchMesh(std::vector<NURBSPatch> &nurbs_patches) {
    patches = nurbs_patches;
}

bool PatchMesh::intersect(Ray &ray, Interaction &interaction) const {
    lbvhHit(interaction, ray);
    return interaction.type != Interaction::Type::NONE;
}

AABB PatchMesh::getPatch(int pos) {
    return patch_AABB[pos];
}

void PatchMesh::genAABB_for_BVH(BVHNode *now) {
    now->aabb = getPatch(now->begin);
    for (int i = 1; i < now->num; ++i) {
        now->aabb = AABB(getPatch(now->begin + i), now->aabb);
    }
}

float PatchMesh::calCost(const AABB& a, int numa, const AABB& b, int numb, const AABB& N) {
    // Calculate the SAH cost
    return a.getSurfaceArea() / N.getSurfaceArea() * numa + b.getSurfaceArea() / N.getSurfaceArea() * numb;
}

int PatchMesh::getPartitionMethod(BVHNode *now) {
    // SAH
    int l = now->begin, r = now->begin + now->num - 1;
    int num = now->num;
    int lsize = 1;
    std::vector<AABB> lprefix_AABB, rprefix_AABB;
    lprefix_AABB.push_back(AABB(getPatch(l)));
    rprefix_AABB.push_back(AABB(getPatch(r)));
    for (int i = 1; i < num; ++i) {
        lprefix_AABB.emplace_back(lprefix_AABB[i - 1], getPatch(l + i));
        rprefix_AABB.emplace_back(rprefix_AABB[i - 1], getPatch(r - i));
    }
    float cost = calCost(lprefix_AABB[0], 1, rprefix_AABB[num - 2], num - 1, now->aabb);
    for (int size = 2; size < num - 1; ++size) {
        float tmpcost = calCost(lprefix_AABB[size - 1], size, rprefix_AABB[num - size - 1], num - size, now->aabb);
        if (cost > tmpcost) {
            cost = tmpcost;
            lsize = size;
        }
    }
    return lsize;
}

void PatchMesh::sort_aabb(int begin, int end, int axis) {
    // Sort the aabb to partition
    if (begin >= end) return;
    // if(debug_flag == 1) printf("%d %d\n", begin, end);
    AABB pivot_aabb = patch_AABB[begin];
    int pivot_index = patch_indices[begin];
    int l = begin, r = end, indx = begin;

    for (; r > l;) {
        for (; r > l;) {
            if (patch_AABB[r].getCenter()[axis] < pivot_aabb.getCenter()[axis]) {
                patch_AABB[l] = patch_AABB[r];
                patch_indices[l] = patch_indices[r];
                indx = r;
                l++;
                break;
            }
            r--;
        }
        for (; r > l;) {
            if (patch_AABB[l].getCenter()[axis] > pivot_aabb.getCenter()[axis]) {
                patch_AABB[r] = patch_AABB[l];
                patch_indices[r] = patch_indices[l];
                indx = l;
                r--;
                break;
            }
            l++;
        }
    }
    patch_AABB[indx] = pivot_aabb;
    patch_indices[indx] = pivot_index;

    sort_aabb(begin, indx - 1, axis);
    sort_aabb(indx + 1, end, axis);
}

void PatchMesh::buildBVH_partition(BVHNode *now, int pre_axis) {
    // TODO partition now node into two child node
    genAABB_for_BVH(now);
    // printf("%d %d\n", now->begin, now->num);
    // printf("%f %f %f\n", now->aabb.getDist(0), now->aabb.getDist(1), now->aabb.getDist(2));
    if (now->num == 1) return;
    float dis = now->aabb.getDist(0);
    now->partition_axis = 0;
    for (int i = 1; i < 3; ++i) {
        if (dis < now->aabb.getDist(i)) {
            dis = now->aabb.getDist(i);
            now->partition_axis = i;
        }
    }

    // sort [begin, begin + num]
    int num = now->num, begin = now->begin;
    // printf("o\n");
    // if(num == 217853 && begin == 0) debug_flag = 1;
    if (pre_axis != now->partition_axis) sort_aabb(begin, begin + num - 1, now->partition_axis);
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
}

int PatchMesh::DFS_BVHTree(BVHNode *now) {
    // TODO : travel the bvh tree to get the DFS order
    lbvh.emplace_back();
    int pos = lbvh.size() - 1;
    if (now->partition_axis == -1) {
        lbvh[pos].set(now, -1);
        return 1;
    }
    int lsize = DFS_BVHTree(now->left);
    int rsize = DFS_BVHTree(now->right);
    lbvh[pos].set(now, lsize + 1);
    return lsize + rsize + 1;
}

void PatchMesh::buildBVH() {
    // TODO ; get the linear bvh
    BVHNode *bvh;
    for (int i = 0; i < patches.size(); ++i) {
        patch_AABB.emplace_back(patches[i]);
        patch_indices.push_back(i);
        auto x = patch_AABB[i];
        // printf("%d low : %f %f %f upper: %f %f %f\n", i, x.low_bnd.x(), x.low_bnd.y(), x.low_bnd.z(),
        // x.upper_bnd.x(), x.upper_bnd.y(), x.upper_bnd.z());
    }
    bvh = new BVHNode();
    bvh->num = (int) patch_AABB.size();
    bvh->begin = 0;
    buildBVH_partition(bvh, -1);
    DFS_BVHTree(bvh);
    delete bvh;
}

void PatchMesh::lbvhHit(Interaction &interaction, Ray &ray) const {
    // TODO: Use the same BVH traversal to traverse the AABB.
    // When reaches a leaf node, call `intersectOnePatch` to determine the nearest intersection.
    std::queue<int> Node_index;
    Node_index.push(0);
    for (; !Node_index.empty();) {
        int now_index = Node_index.front();
        Node_index.pop();
        auto now = lbvh[now_index];
        float t_in, t_out;
        if (now.aabb.intersect(ray, &t_in, &t_out)) {
            if (now.num != -1) {
                for (int j = 0; j < now.num; ++j) {
                    int i = patch_indices[now.object_offset + j];
                    bool tag = false;
                    auto patch = patches[i];
                    int step = 10;
                    Interaction temp1, temp2;
                    float initial_u = (patches[i].range_u_.x() + patches[i].range_u_.y()) / 2.0f;
                    float initial_v = (patches[i].range_v_.x() + patches[i].range_v_.y()) / 2.0f;
                    if (intersectOnePatch(ray, temp2, patches[i], Vec2f(initial_u, initial_v)) &&
                        temp2.dist < interaction.dist)
                        interaction = temp2;
                    // if (intersectOnePatch(ray, temp, patches[i], Vec2f(0.1f, 0.5f))) {
                    //   // printf("[%f %f] %f %f\n", t_in, t_out, temp.dist, interaction.dist);
                    //   if (temp.dist < interaction.dist) interaction = temp;
                    // }
                }
            } else {
                Node_index.push(now_index + 1);
                Node_index.push(now_index + now.node_offset);
            }
        }
    }
    }

bool PatchMesh::intersectOnePatch(Ray &ray, Interaction &interaction, const NURBSPatch &patch,
                                  Vec2f initial_val) const {
    // Write the ray (o+td) as an intersection of two planes.
    Vec3f N1, N2;
    if (std::abs(ray.direction.x()) > std::abs(ray.direction.y()) &&
        std::abs(ray.direction.x()) > std::abs(ray.direction.z())) {
        N1 = Vec3f(ray.direction.y(), -ray.direction.x(), 0);
    } else {
        N1 = Vec3f(0, ray.direction.z(), -ray.direction.y());
    }
    N1.normalize();
    N2 = N1.cross(ray.direction).normalized();

    float d1 = -N1.dot(ray.origin);
    float d2 = -N2.dot(ray.origin);

    constexpr float eps = EPS;
    // constexpr float eps = 0.001;
    constexpr int MAX_ITER = 7;

    // TODO: write Newton Root-Finder here?
    float error_prev = std::numeric_limits<float>::max();
    // u, v are the initial guess of the intersection point
    // const float u_initial = (patch.range_u_.x() + patch.range_u_.y()) / 2.0f, v_initial = (patch.range_v_.x() +
    // patch.range_v_.y())/2.0f;
    const float u_initial = initial_val.x(), v_initial = initial_val.y();
    float u = u_initial, v = v_initial;
    // printf("(u, v) : %f %f\n", u, v);
    for (int iter = 0; iter < MAX_ITER; ++iter) {
        // S = evaluate surface (u, v)
        // printf("(u, v) : %f %f\n", u, v);

        auto res = patch.evaluate(u, v);
        Vertex S = res.first;
        // printf("(%f, %f, %f) u: %f v : %f\n", S.position.x(), S.position.y(), S.position.z(), u, v);
        Vec2f R = Vec2f(N1.dot(S.position) + d1, N2.dot(S.position) + d2);
        float error = R.norm();
        // printf("%f\n", error);
        if (error < eps) {
            // printf("in\n");

            // Compute the intersection point
            // Update the interaction variables
            float t = (S.position - ray.origin).dot(ray.direction);
            if (t < EPS) return false;
            interaction.pos = S.position;
            interaction.dist = t;
            // printf("S position: %f %f %f\n", S.position.x(), S.position.y(), S.position.z());
            // printf("o: %f %f %f\n", ray.origin.x(), ray.origin.y(), ray.origin.z());
            // interaction.normal = (patch.evaluate(u+eps, v).position - patch.evaluate(u-eps,
            // v).position).cross(patch.evaluate(u, v+eps).position - patch.evaluate(u, v-eps).position).normalized();
            interaction.normal = S.normal;
            // printf("(%f %f) pos: %f %f %f\n",u, v, S.position.x(), S.position.y(), S.position.z());
            // auto ref = Vec3f(S.position.x() + 0.3f, S.position.y() - 1.0f, S.position.z() - 0.4f).normalized();
            // printf("pos - o: %f %f %f\n", ref.x(), ref.y(), ref.z());
            // printf("normal: %f %f %f\n", S.normal.x(), S.normal.y(), S.normal.z());
            interaction.material = bsdf;
            // interaction.wi = ray.direction;
            interaction.wo = -ray.direction;
            interaction.type = Interaction::Type::NURBS;
            // printf("Patch mesh : %f\n", interaction.dist);
            return true;
        }
        // printf("%f\n", error);
        if (std::abs(error) > error_prev) {
            // printf("\n");
            // printf("out (error increase)\n");
            return false;
        }
        if (u < patch.range_u_.x() || u >= patch.range_u_.y()) {
            // printf("out (out of range u)\n");
            return false;
        }
        if (v < patch.range_v_.x() || v >= patch.range_v_.y()) {
            // printf("out (out of range v)\n");
            return false;
        }
        error_prev = std::abs(error);
        // J = compute Jacobian matrix
        auto ress = patch.evaluate(u, v).second;
        Vec2f Fu = {N1.dot(ress.first), N2.dot(ress.first)};
        Vec2f Fv = {N1.dot(ress.second), N2.dot(ress.second)};
        Mat2f J;
        // printf("%f %f\n", Fu.x(), Fu.y());
        // printf("%f %f\n", Fv.x(), Fv.y());
        J.col(0) = Fu;
        J.col(1) = Fv;
        // printf("%f %f %f %f\n", J(0, 0), J(1, 0), J(1, 0), J(0, 0));

        // if J is singular
        if (std::abs(J.determinant()) < eps) {
            u += 0.1 * (u_initial - u) * drand48();
            v += 0.1 * (v_initial - v) * drand48();
        } else {
            // u -= J.inverse().transpose().col(0).dot(R);
            // v -= J.inverse().transpose().col(1).dot(R);
            u -= J.inverse().row(0).dot(R);
            v -= J.inverse().row(1).dot(R);
        }
    }
    // printf("out (out of itr time)\n");
    // printf("\n");
    return false;
}

