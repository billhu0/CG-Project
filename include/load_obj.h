#ifndef LOAD_OBJ_H_
#define LOAD_OBJ_H_

#include "core.h"
#include "geometry.h"
#include "nurbs.h"
#include <string>

static bool loadObj(const std::string &path, std::vector<Vec3f> &vertices, std::vector<Vec3f> &normals,
                    std::vector<int> &v_index, std::vector<int> &n_index);

std::shared_ptr<TriangleMesh> makeMeshObject(std::string path_to_obj, Vec3f translation, float scale);

std::vector<NURBSSurface> loadNURBS(const std::string &path_to_obj);

std::vector<std::shared_ptr<Mesh>> makeNURBSObject(const std::string &path_to_obj, const Vec3f &translation,
                                                   float scale);

#endif  // LOAD_OBJ_H_
