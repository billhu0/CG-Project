#include "nurbs.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>


// bool 
// samen(Triangle x, Triangle y)
// {
//     return ((x.p1 == y.p1) && (x.p2 == y.p2) && (x.p3 == y.p3));
// }

int 
NURBSCurve::Findk(float t)
{
    for(int i = 0; i < knots_.size(); ++ i)
    {
        if(t == 1.0f && knots_[i + 1] == 1.0f) return i;
        if(knots_[i] <= t && knots_[i + 1] > t) return i;
    }
}

NURBSCurve::NURBSCurve(std::vector<Vec3f> control_points, std::vector<float> w, std::vector<float> knots, int degree)
{
    control_points_ = control_points;
    w_ = w;
    knots_ = knots;
    degree_ = degree;
}

NURBSCurve::NURBSCurve(int n, int degree)
{
    control_points_.resize(n);
    w_.resize(n);
    degree_ = degree;
    knots_.resize(n + degree + 1);
}

void 
NURBSCurve::setControlPointAndWeight(int i, Vec3f control_ponit, float w)
{
    control_points_[i] = control_ponit;
    w_[i] = w;
}

void
NURBSCurve::setKnot(int i, float knot)
{
    knots_[i] = knot;
}

void
NURBSCurve::setKnots(std::vector<float> knots)
{
    knots_ = knots;
}

void
NURBSCurve::setKnotUniform()
{
    int n = control_points_.size();
    float delta = 1.0f / (float)(n - degree_);
    for(int i = 0; i < knots_.size(); ++ i)
    {
        if(i < degree_) knots_[i] = (float)0.0;
        else if(i > n) knots_[i] = (float)1.0;
        else knots_[i] = (float)(i - degree_) * delta;
    }
}

std::pair<Vertex, float>
NURBSCurve::evaluate(float t)
{
    int k = Findk(t);
    std::vector<Vec4f> P;
    Vertex res;
    float w;
    for(int i = 0; i < degree_ + 1; ++ i)
    {
        P.push_back(Vec4f(control_points_[i + k - degree_].x(), control_points_[i + k - degree_].y(), control_points_[i + k - degree_].z(), 1.0f) * w_[i + k - degree_]);
    }
    for(int j = 1; j <= degree_; ++ j)
    {
        for(int i = degree_; i >= j; -- i)
        {
            float alpha = (t - knots_[k + i - degree_]) / (knots_[k + i - (j - 1)] -  knots_[k + i - degree_]);
            P[i] = (1.0f - alpha) * P[i - 1] + alpha * P[i];
        }
        if(j == degree_ - 1) res.normal = Vec3f(P[degree_].x() / P[degree_].w(), P[degree_].y() / P[degree_].w(), P[degree_].z() / P[degree_].w()) - Vec3f(P[degree_ - 1].x() / P[degree_ - 1].w(), P[degree_ - 1].y() / P[degree_ - 1].w(), P[degree_ - 1].z() / P[degree_ - 1].w());
        if(j == degree_)
        {
            res.position = Vec3f(P[degree_].x() / P[degree_].w(), P[degree_].y() / P[degree_].w(), P[degree_].z() / P[degree_].w());
            w = P[degree_].w();
        }
    }
    return std::pair(res, w);
}

NURBSSurface::NURBSSurface(int m, int n, int degreem, int degreen)
{
    control_points_m_.resize(m);
    w_m_.resize(m);
    for(auto &sub_vec : control_points_m_) sub_vec.resize(n);
    for(auto &sub_vec : w_m_) sub_vec.resize(n);
    degree_m_ = degreem;
    knots_m_.resize(m + degreem + 1);
    control_points_n_.resize(n);
    w_n_.resize(n);
    for(auto &sub_vec : control_points_n_) sub_vec.resize(m);
    for(auto &sub_vec : w_n_) sub_vec.resize(m);
    degree_n_ = degreen;
    knots_n_.resize(n + degreen + 1);
}

void 
NURBSSurface::setControlPointAndWeight(int i, int j, Vec3f control_point, float w)
{
    control_points_m_[i][j] = control_point;
    control_points_n_[j][i] = control_point;
    w_m_[i][j] = w;
    w_n_[j][i] = w;
}

void
NURBSSurface::setKnotM(int i, float knot)
{
    knots_m_[i] = knot;
}

void
NURBSSurface::setKnotN(int i, float knot)
{
    knots_n_[i] = knot;
}

Vertex
NURBSSurface::evaluate(float u, float v)
{
    Vertex point_m, point_n, point;
    NURBSCurve curve_m(control_points_m_.size(), degree_m_);
    curve_m.setKnots(knots_m_);
    for(int i = 0; i < control_points_m_.size(); ++ i)
    {
        NURBSCurve tmp(control_points_m_[i], w_m_[i], knots_n_, degree_n_);
        auto tmp_retval = tmp.evaluate(v);
        curve_m.setControlPointAndWeight(i, tmp_retval.first.position, tmp_retval.second);
    }
    point_m = curve_m.evaluate(u).first;
    NURBSCurve curve_n(control_points_n_.size(), degree_n_);
    curve_n.setKnots(knots_n_);
    for(int i = 0; i < control_points_n_.size(); ++ i)
    {
        NURBSCurve tmp(control_points_n_[i], w_n_[i], knots_m_, degree_m_);
        auto tmp_retval = tmp.evaluate(u);
        curve_n.setControlPointAndWeight(i, tmp_retval.first.position, tmp_retval.second);
    }
    point_n = curve_n.evaluate(v).first;
    point.position = (point_m.position + point_n.position) * 0.5f;
    point.normal = point_m.normal.cross(point_n.normal).normalized();
    return point;
}

std::shared_ptr<TriangleMesh> 
NURBSSurface::genObject_uniformsample(const Vec3f& translation, float scale)
{
    std::vector<Vec3f> vertices;
    std::vector<Vec3f> normals;
    std::vector<int> v_idx;
    std::vector<int> n_idx;

    unsigned int divide_num = 100;
    unsigned int line_size = divide_num;
    
    for(int i = 0; i <= divide_num; ++ i)
    {
        float u = (float) i / (float) divide_num;
        
        for(int j = 0; j <= divide_num; ++ j)
        {
            float v = (float) j / (float) divide_num;
            auto res = evaluate(u, v);
            vertices.push_back(res.position * scale + translation);
            // vertices.push_back(res.position);
            // if(res.normal.norm() < 0.01f) res.normal = Vec3f(0.000000f, 0.219522f, -0.975608f);
            normals.push_back(res.normal);
            // std::cout << evaluate(u, v).position.x << " " << evaluate(u, v).position.y << " " << evaluate(u, v).position.z << "\n";

        }
    }
    for(unsigned int i = 0; i < line_size; ++ i)
        for(unsigned int j = 0; j < line_size; ++ j)
        {
            v_idx.push_back(i * (line_size + 1) + j);
            v_idx.push_back(i * (line_size + 1) + j + 1);
            v_idx.push_back((i + 1) * (line_size + 1) + j);
            v_idx.push_back(i * (line_size + 1) + j + 1);
            v_idx.push_back((i + 1) * (line_size + 1) + j);
            v_idx.push_back((i + 1) * (line_size + 1) + j + 1);

            n_idx.push_back(i * (line_size + 1) + j);
            n_idx.push_back(i * (line_size + 1) + j + 1);
            n_idx.push_back((i + 1) * (line_size + 1) + j);
            n_idx.push_back(i * (line_size + 1) + j + 1);
            n_idx.push_back((i + 1) * (line_size + 1) + j);
            n_idx.push_back((i + 1) * (line_size + 1) + j + 1);
        }
    // for (auto v : vertices) printf("v %f %f %f\n", v.x(), v.y(), v.z());
    // for (auto n : normals) printf("vn %f %f %f\n", n.x(), n.y(), n.z());
    // printf ("vt 0 1\n");
    // int num = 0;
    // // printf("v_indx : ");
    // for (auto index : v_idx)
    // {
    //     if(num == 0) printf("f ");
    //     printf("%d/1/%d ", index + 1, index + 1);
    //     num ++;
    //     if(num == 3)
    //     {
    //         printf("\n");
    //         num = 0;
    //     }
    // }
//   for (auto v : vertices) printf("v %f %f %f\n", v.x(), v.y(), v.z());
//   for (auto n : normals) printf("n %f %f %f\n", n.x(), n.y(), n.z());
//   printf("v_indx : ");
//   for (auto index : v_idx) printf ("%d ", index);
//   printf("\nn_indx :");
//   for (auto index : n_idx) printf ("%d ", index);
    auto ptr = std::make_shared<TriangleMesh>(vertices, normals, v_idx, n_idx);
    // ptr->print_triangle_mesh();
    return ptr;
}



void
NURBSSurface::setKnotMUniform()
{
    int n = control_points_m_.size();
    float delta = 1.0f / (float)(n - degree_m_);
    for(int i = 0; i < knots_m_.size(); ++ i)
    {
        if(i < degree_m_) knots_m_[i] = (float)0.0;
        else if(i > n) knots_m_[i] = (float)1.0;
        else knots_m_[i] = (float)(i - degree_m_) * delta;
    }
}

void
NURBSSurface::setKnotNUniform()
{
    int n = control_points_n_.size();
    float delta = 1.0f / (float)(n - degree_n_);
    for(int i = 0; i < knots_n_.size(); ++ i)
    {
        if(i < degree_n_) knots_n_[i] = (float)0.0;
        else if(i > n) knots_n_[i] = (float)1.0;
        else knots_n_[i] = (float)(i - degree_n_) * delta;
    }
}

// Object 
// NURBSSurface::genObject_adaptivesample(float tolerance)
// {
//     Object res;
//     std::vector<std::pair<float, float>> parameter;
//     int triangle_num = 0;
//     unsigned int divide_num = 50;
//     unsigned int vertex_num = 0;
//     unsigned int line_size = divide_num;
//     for(int i = 0; i <= divide_num; ++ i)
//     {
//         float u = (float) i / (float) divide_num;
        
//         for(int j = 0; j <= divide_num; ++ j)
//         {
//             float v = (float) j / (float) divide_num;
//             res.vertices.push_back(evaluate(u, v));
//             parameter.push_back(std::make_pair(u, v));
//             vertex_num ++;
//         }
//     }
//     std::vector<Triangle> T;
//     for(unsigned int i = 0; i < line_size; ++ i)
//         for(unsigned int j = 0; j < line_size; ++ j)
//         {
//             T.push_back(Triangle(i * (line_size + 1) + j, i * (line_size + 1) + j + 1, (i + 1) * (line_size + 1) + j));
//             T.push_back(Triangle(i * (line_size + 1) + j + 1, (i + 1) * (line_size + 1) + j, (i + 1) * (line_size + 1) + j + 1));
//         }
//     std::vector<Triangle> NT;
//     for(;T.size() != NT.size();) //there is no subdivide in the last loop
//     {
//         if(!NT.empty())
//         {
//             T.clear();
//             T = NT;
//             NT.clear();
//         }
//         for(auto t : T)
//         {
//             Vertex p1 = res.vertices[t.p1], p2 = res.vertices[t.p2], p3 = res.vertices[t.p3];
//             std::pair parameter1 = parameter[t.p1], parameter2 = parameter[t.p2], parameter3 = parameter[t.p3];
//             if((!t.flag12) && (!check_flat(p1.normal, p2.normal, tolerance)))
//             {
//                 t.flag12 = true;
//                 t.div_edge_num ++;
//                 float u = (parameter1.first + parameter2.first) * 0.5f;
//                 float v = (parameter1.second + parameter2.second) * 0.5f;
//                 res.vertices.push_back(evaluate(u, v));
//                 parameter.push_back(std::make_pair(u, v));
//                 t.p12 = vertex_num;
//                 for(auto tt : T)
//                 {
//                     if(samen(tt, t)) continue;
//                     if((tt.p1 == t.p1 && tt.p2 == t.p2)||(tt.p1 == t.p2 && tt.p2 == t.p1)) {tt.flag12 = true; tt.div_edge_num++; tt.p12 = vertex_num; break;}
//                     if((tt.p2 == t.p1 && tt.p3 == t.p2)||(tt.p2 == t.p2 && tt.p3 == t.p1)) {tt.flag23 = true; tt.div_edge_num++; tt.p23 = vertex_num; break;}
//                     if((tt.p3 == t.p1 && tt.p1 == t.p2)||(tt.p3 == t.p2 && tt.p1 == t.p1)) {tt.flag31 = true; tt.div_edge_num++; tt.p31 = vertex_num; break;}
//                 }
//                 vertex_num ++;
//             }
//             if((!t.flag23) && (!check_flat(p2.normal, p3.normal, tolerance)))
//             {
//                 t.flag23 = true;
//                 t.div_edge_num ++;
//                 float u = (parameter2.first + parameter3.first) * 0.5f;
//                 float v = (parameter2.second + parameter3.second) * 0.5f;
//                 res.vertices.push_back(evaluate(u, v));
//                 parameter.push_back(std::make_pair(u, v));
//                 t.p23 = vertex_num;
//                 for(auto tt : T)
//                 {
//                     if(samen(tt, t)) continue;
//                     if((tt.p1 == t.p2 && tt.p2 == t.p3)||(tt.p1 == t.p3 && tt.p2 == t.p2)) {tt.flag12 = true; tt.div_edge_num++; tt.p12 = vertex_num; break;}
//                     if((tt.p2 == t.p2 && tt.p3 == t.p3)||(tt.p2 == t.p3 && tt.p3 == t.p2)) {tt.flag23 = true; tt.div_edge_num++; tt.p23 = vertex_num; break;}
//                     if((tt.p3 == t.p2 && tt.p1 == t.p3)||(tt.p3 == t.p3 && tt.p1 == t.p2)) {tt.flag31 = true; tt.div_edge_num++; tt.p31 = vertex_num; break;}
//                 }
//                 vertex_num ++;
//             }
//             if((!t.flag31) && (!check_flat(p3.normal, p1.normal, tolerance)))
//             {
//                 t.flag31 = true;
//                 t.div_edge_num ++;
//                 float u = (parameter3.first + parameter1.first) * 0.5f;
//                 float v = (parameter3.second + parameter1.second) * 0.5f;
//                 res.vertices.push_back(evaluate(u, v));
//                 parameter.push_back(std::make_pair(u, v));
//                 t.p31 = vertex_num;
//                 for(auto tt : T)
//                 {
//                     if(samen(tt, t)) continue;
//                     if((tt.p1 == t.p3 && tt.p2 == t.p1)||(tt.p1 == t.p1 && tt.p2 == t.p3)) {tt.flag12 = true; tt.div_edge_num++; tt.p12 = vertex_num; break;}
//                     if((tt.p2 == t.p3 && tt.p3 == t.p1)||(tt.p2 == t.p1 && tt.p3 == t.p3)) {tt.flag23 = true; tt.div_edge_num++; tt.p23 = vertex_num; break;}
//                     if((tt.p3 == t.p3 && tt.p1 == t.p1)||(tt.p3 == t.p1 && tt.p1 == t.p3)) {tt.flag31 = true; tt.div_edge_num++; tt.p31 = vertex_num; break;}
//                 }
//                 vertex_num ++;
//             }

//             if(t.div_edge_num == 0)
//             {
//                 NT.push_back(t);
//             }
//             else if(t.div_edge_num == 1)
//             {
//                 if(t.flag12)
//                 {
//                     NT.push_back(Triangle(t.p1, t.p3, t.p12));
//                     NT.push_back(Triangle(t.p12, t.p2, t.p3));
//                 }
//                 else if(t.flag23)
//                 {
//                     NT.push_back(Triangle(t.p1, t.p3, t.p23));
//                     NT.push_back(Triangle(t.p1, t.p2, t.p23)); 
//                 }
//                 else
//                 {
//                     NT.push_back(Triangle(t.p2, t.p3, t.p31));
//                     NT.push_back(Triangle(t.p1, t.p2, t.p31));   
//                 }
//             }
//             else if(t.div_edge_num == 2)
//             {
//                 if(!t.flag12)
//                 {
//                     NT.push_back(Triangle(t.p3, t.p31, t.p23));
//                     NT.push_back(Triangle(t.p1, t.p31, t.p23));
//                     NT.push_back(Triangle(t.p1, t.p2, t.p23));
//                 }
//                 else if(!t.flag23)
//                 {
//                     NT.push_back(Triangle(t.p2, t.p31, t.p12));
//                     NT.push_back(Triangle(t.p1, t.p31, t.p12));
//                     NT.push_back(Triangle(t.p3, t.p2, t.p31)); 
//                 }
//                 else
//                 {
//                     NT.push_back(Triangle(t.p2, t.p23, t.p12));
//                     NT.push_back(Triangle(t.p3, t.p23, t.p12));
//                     NT.push_back(Triangle(t.p3, t.p1, t.p12));  
//                 }
//             }
//             else
//             {
//                 NT.push_back(Triangle(t.p3, t.p31, t.p23));
//                 NT.push_back(Triangle(t.p1, t.p31, t.p12));
//                 NT.push_back(Triangle(t.p31, t.p23, t.p12));
//                 NT.push_back(Triangle(t.p2, t.p23, t.p12));
//             }
//         }
//     }   

//     for(auto t : NT)
//     {
//         res.indices.push_back(t.p1);
//         res.indices.push_back(t.p2);
//         res.indices.push_back(t.p3);
//     }
//     res.draw_mode.drawmethod = DRAW_ELEMENTS_WITH_SHADER;
//     res.draw_mode.primitive_mode = GL_TRIANGLES;


//     res.init();
//     return res;
// }

// bool
// NURBSSurface::check_flat(Vec3f normal1, Vec3f normal2, float tolerance)
// {
//     return (1 - glm::dot(normal1, normal2)) < tolerance;
// }
