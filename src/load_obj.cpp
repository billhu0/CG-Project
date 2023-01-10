#include "load_obj.h"

#define TINYOBJLOADER_IMPLEMENTATION
#include <tiny_obj_loader.h>
#include <iostream>
#include "nurbs.h"

static bool loadObj(const std::string &path, std::vector<Vec3f> &vertices, std::vector<Vec3f> &normals,
                    std::vector<int> &v_index, std::vector<int> &n_index) {
    std::cout << "-- Loading model " << path << std::endl;

    tinyobj::ObjReaderConfig readerConfig;
    // readerConfig.mtl_search_path = "./";  // Path to material files

    tinyobj::ObjReader reader;

    if (!reader.ParseFromFile(path, readerConfig)) {
        if (!reader.Error().empty()) {
            std::cerr << "TinyObjReader: " << reader.Error();
        }
        exit(1);
    }

    if (!reader.Warning().empty()) {
        std::cout << "TinyObjReader: " << reader.Warning();
    }

    auto &attrib = reader.GetAttrib();
    auto &shapes = reader.GetShapes();
    auto &materials = reader.GetMaterials();

    for (size_t i = 0; i < attrib.vertices.size(); i += 3) {
        vertices.emplace_back(attrib.vertices[i], attrib.vertices[i + 1], attrib.vertices[i + 2]);
    }
    // for (size_t i = 0; i < attrib.texcoords.size(); i += 2) {
    //   texCoords.push_back(vec2(attrib.texcoords[i], attrib.texcoords[i + 1]));
    // }
    for (size_t i = 0; i < attrib.normals.size(); i += 3) {
        normals.emplace_back(attrib.normals[i], attrib.normals[i + 1], attrib.normals[i + 2]);
    }
    // Loop over shapes
    for (size_t s = 0; s < shapes.size(); s++) {
        // Loop over faces(polygon)
        size_t index_offset = 0;
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
            auto fv = size_t(shapes[s].mesh.num_face_vertices[f]);

            // Loop over vertices in the face.
            for (size_t v = 0; v < fv; v++) {
                tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
                v_index.push_back(idx.vertex_index);
                // tIndex.push_back(idx.texcoord_index);
                n_index.push_back(idx.normal_index);
            }
            index_offset += fv;
        }
    }
    // std::cout << "  # vertices: " << attrib.vertices.size() / 3 << std::endl;
    // std::cout << "  # faces: " << v_index.size() / 3 << std::endl;
    return true;
}

std::shared_ptr<TriangleMesh> makeMeshObject(std::string path_to_obj, Vec3f translation, float scale) {
    std::vector<Vec3f> vertices;
    std::vector<Vec3f> normals;
    std::vector<int> v_idx;
    std::vector<int> n_idx;
    loadObj(path_to_obj, vertices, normals, v_idx, n_idx);
    for (auto &v : vertices) v = v * scale + translation;
    // printf("\n%s\n", path_to_obj.c_str());
    // for (auto v : vertices) printf("v %f %f %f\n", v.x(), v.y(), v.z());
    // for (auto n : normals) printf("n %f %f %f\n", n.x(), n.y(), n.z());
    // printf("v_indx : ");
    // for (auto index : v_idx) printf ("%d ", index);
    // printf("\nn_indx :");
    // for (auto index : n_idx) printf ("%d ", index);
    auto ptr = std::make_shared<TriangleMesh>(vertices, normals, v_idx, n_idx);
    // ptr->print_triangle_mesh();
    return ptr;
}

std::vector<NURBSSurface> loadNURBS(const std::string &path_to_obj) {
    std::ifstream infile(path_to_obj);
    // Check the file open successfully
    if (!infile.is_open()) {
        std::cout << "Error opening the file " << path_to_obj << std::endl;
        exit(1);
    }
    unsigned int surface_num;
    unsigned int m, n, dgree_m, dgree_n;
    std::vector<std::vector<unsigned int>> indecies;
    std::vector<Vec3f> control_point_pos;
    std::string line;
    std::getline(infile, line);
    std::istringstream istr_all(line);
    istr_all >> surface_num;
    std::vector<NURBSSurface> retval;
    for (int i = 0; i < surface_num; ++i) {
        std::getline(infile, line);
        std::istringstream istr(line);
        istr >> m >> n >> dgree_m >> dgree_n;
        NURBSSurface surface(m, n, dgree_m, dgree_n);
        std::getline(infile, line);
        istr = std::istringstream(line);
        for (int j = 0; j < m + dgree_m + 1; ++j) {
            float knot;
            istr >> knot;
            surface.setKnotM(j, knot);
        }
        std::getline(infile, line);
        istr = std::istringstream(line);
        for (int j = 0; j < n + dgree_n + 1; ++j) {
            float knot;
            istr >> knot;
            surface.setKnotN(j, knot);
        }
        for (int j = 0; j < m * n; ++j) {
            std::getline(infile, line);
            istr = std::istringstream(line);
            int indxi, indxj;
            float x, y, z, w;
            istr >> indxi >> indxj >> x >> y >> z >> w;
            surface.setControlPointAndWeight(indxi, indxj, Vec3f(x, y, z), w);
        }
        retval.push_back(surface);
    }
    return retval;
}

std::vector<std::shared_ptr<Mesh>> makeNURBSObject(const std::string &path_to_obj, const Vec3f &translation,
                                                   float scale) {
    std::vector<std::shared_ptr<Mesh>> retval;
    auto NURBS_objs = loadNURBS(path_to_obj);
    for (auto nurbs : NURBS_objs) {
        if (ENABLE_NURBS)
            retval.push_back(nurbs.genMesh_triangle(translation, scale));
        else
            retval.push_back(nurbs.genMesh_patch(translation, scale));
    }
    return retval;
}
