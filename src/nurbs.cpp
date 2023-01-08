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
NURBSSurface::genMesh_triangle(const Vec3f& translation, float scale)
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

std::shared_ptr<PatchMesh> 
NURBSSurface::genMesh_patch(const Vec3f& translation, float scale)
{

    float C = 1.0f;
    //generate V and A
    std::vector<std::vector<Vec3f>> V_m, V_n, A_m, A_n;
    int m = control_points_m_.size(), n = control_points_n_.size();
    V_m.resize(m);  A_m.resize(m);
    V_n.resize(n);  A_n.resize(n);
    for(int i = 0; i < m; ++ i)
    {
        V_m[i].resize(n);
        V_m[i][0] = Vec3f(0, 0, 0);
        A_m[i].resize(n);
        A_m[i][0] = Vec3f(0, 0, 0);
        for(int j = 1; j < n; ++ j)
        {
            V_m[i][j] = degree_n_ * (control_points_m_[i][j] - control_points_m_[i][j - 1]) / (knots_n_[j + degree_n_] - knots_n_[j]);
            A_m[i][j] = (degree_n_ - 1) * (V_m[i][j] - V_m[i][j - 1]) / (knots_n_[j + degree_n_ - 1] - knots_n_[j]);
        }
    }
    for(int i = 0; i < n; ++ i)
    {
        V_n[i].resize(m);
        V_n[i][0] = Vec3f(0, 0, 0);
        A_n[i].resize(m);
        A_n[i][0] = Vec3f(0, 0, 0);
        for(int j = 1; j < m; ++ j)
        {
            V_n[i][j] = degree_m_ * (control_points_n_[i][j] - control_points_n_[i][j - 1]) / (knots_m_[j + degree_m_] - knots_m_[j]);
            A_n[i][j] = (degree_m_ - 1) * (V_n[i][j] - V_n[i][j - 1]) / (knots_m_[j + degree_m_ - 1] - knots_m_[j]);
        }
    }

    //caculate the add num
    std::vector<int> knots_m_addnum, knots_n_addnum;
    int m_add_sum = 0, n_add_sum = 0;
    knots_m_addnum.resize(knots_m_.size() - 1);
    knots_n_addnum.resize(knots_n_.size() - 1);
    for(int i = 0; i < knots_m_.size() - 1; ++ i)
    {
        float delta = knots_m_[i + 1] - knots_m_[i];
        knots_m_addnum[i] = 0;
        if(delta < EPS) continue;
        for(int j = 0; j < knots_n_.size() - 1; ++ j)
        {
            if(knots_n_[j + 1] == knots_n_[j]) continue;
            float maxa = 0, sumv = V_n[j][i - degree_m_ + 1].norm();
            for(int k = i - degree_m_ + 2; k <= i; ++ k)
            {
                auto tmp = A_n[j][k].norm();
                maxa = maxa > tmp ? maxa : tmp;
                sumv += V_n[j][k].norm();
            }
            int N = C * maxa * powf(delta, 1.5f) / powf(sumv / degree_m_, 0.5f);
            knots_m_addnum[i] = knots_m_addnum[i] > N ? knots_m_addnum[i] : N;
        }
        m_add_sum = m_add_sum + knots_m_addnum[i];
    }
    for(int i = 0; i < knots_n_.size() - 1; ++ i)
    {
        float delta = knots_n_[i + 1] - knots_n_[i];
        knots_n_addnum[i] = 0;
        if(delta < EPS) continue;
        for(int j = 0; j < knots_m_.size() - 1; ++ j)
        {
            if(knots_m_[j + 1] == knots_m_[j]) continue;
            float maxa = 0, sumv = V_m[j][i - degree_n_ + 1].norm();
            for(int k = i - degree_n_ + 2; k <= i; ++ k)
            {
                auto tmp = A_m[j][k].norm();
                maxa = maxa > tmp ? maxa : tmp;
                sumv = V_m[j][k].norm();
            }
            int N = C * maxa * powf(delta, 1.5f) / powf(sumv / degree_n_, 0.5f);
            knots_n_addnum[i] = knots_n_addnum[i] > N ? knots_n_addnum[i] : N;
        }
        n_add_sum = n_add_sum + knots_n_addnum[i];
    }

    //subdivision
    //refine

    //refine the n direction
    std::vector<std::vector<Vec3f>> new_control_points_n = control_points_n_;
    std::vector<std::vector<float>> new_w_n = w_n_;
    std::vector<float> new_knots_m_;
    new_knots_m_.resize(knots_m_.size() + m_add_sum * degree_m_);
    int new_knots_length_m = 0;
    for(int i = 0; i < knots_m_.size(); ++ i) new_knots_m_[i] = knots_m_[i];
    new_knots_length_m = knots_m_.size();
    for(int ii = 0; ii < knots_m_.size() - 1; ++ ii) 
    {
        if(knots_m_[ii] == knots_m_[ii + 1]) continue;
        int vd = knots_m_addnum[ii] + 1;
        float delta = (knots_m_[ii + 1] - knots_m_[ii]) / vd;
        for(int j = 1; j <= knots_m_addnum[ii]; ++ j)
        {
            float t = knots_m_[ii] + j * delta;
            int i = 0;
            for(;i + 1 < new_knots_length_m && (!(new_knots_m_[i] <= t && new_knots_m_[i + 1] > t)); ++ i);
            for(int k = 0; k < control_points_n_.size(); ++ k)
            {
                std::vector<Vec3f> tmp_controlpoints;
                auto sizet = new_control_points_n[k].size() + 1;
                tmp_controlpoints.resize(sizet);
                std::vector<float> tmp_weight;
                tmp_weight.resize(sizet);

                for(int l = 0; l <= i - degree_m_; ++ l)
                {
                    tmp_controlpoints[l] = new_control_points_n[k][l];
                    tmp_weight[l] = new_w_n[k][l];
                }
                for(int l = i - degree_m_ + 1; l <= i; ++ l)
                {
                    auto alpha = (t - new_knots_m_[l]) / (new_knots_m_[l + degree_m_] - new_knots_m_[l]);
                    tmp_controlpoints[l] = (1.0f - alpha) * new_control_points_n[k][l - 1] + alpha * new_control_points_n[k][l];
                    tmp_weight[l] = (1.0f - alpha) * new_w_n[k][l - 1] + alpha * new_w_n[k][l];
                }
                for(int l = i + 1; l < sizet; ++ l)
                {
                    tmp_controlpoints[l] = new_control_points_n[k][l - 1];
                    tmp_weight[l] = new_w_n[k][l - 1];                    
                }

                new_control_points_n[k] = tmp_controlpoints;
                new_w_n[k] = tmp_weight;
            }
            for(int k = new_knots_length_m - 1; k >= 0; -- k)
            {
                if(new_knots_m_[k] > t) new_knots_m_[k + 1] = new_knots_m_[k];
                else
                {
                    new_knots_m_[k + 1] = t;
                    break;
                }
            }
            new_knots_length_m ++;
        }
    }
    // for(auto x : new_knots_m_) printf("%f ", x);
    // printf("\n");
    //close n direction
    int length_m = new_knots_length_m;
    auto x_knots_m = new_knots_m_;
    for(int ii = degree_m_ + 1; ii < length_m - 1 - degree_m_; ++ ii)
    {
        for(int j = 0; j < degree_m_ - 1; ++ j)
        {
            float t = x_knots_m[ii];
            // printf("%f\n", t);
            for(int k = 0; k < control_points_n_.size(); ++ k)
            {
                std::vector<Vec3f> tmp_controlpoints;
                auto sizet = new_control_points_n[k].size() + 1;
                tmp_controlpoints.resize(sizet);
                std::vector<float> tmp_weight;
                tmp_weight.resize(sizet);
                int i = 0;
                for(;i + 1 < new_knots_length_m && (!(new_knots_m_[i] <= t && new_knots_m_[i + 1] > t)); ++ i);
                // printf("%d\n",i);
                for(int l = 0; l <= i - degree_m_; ++ l)
                {
                    tmp_controlpoints[l] = new_control_points_n[k][l];
                    tmp_weight[l] = new_w_n[k][l];
                }
                for(int l = i - degree_m_ + 1; l <= i; ++ l)
                {
                    auto alpha = (t - new_knots_m_[l]) / (new_knots_m_[l + degree_m_] - new_knots_m_[l]);
                    tmp_controlpoints[l] = (1.0f - alpha) * new_control_points_n[k][l - 1] + alpha * new_control_points_n[k][l];
                    tmp_weight[l] = (1.0f - alpha) * new_w_n[k][l - 1] + alpha * new_w_n[k][l];
                }
                for(int l = i + 1; l < sizet; ++ l)
                {
                    tmp_controlpoints[l] = new_control_points_n[k][l - 1];
                    tmp_weight[l] = new_w_n[k][l - 1];                    
                }

                new_control_points_n[k] = tmp_controlpoints;
                new_w_n[k] = tmp_weight;
            }
            for(int k = new_knots_length_m - 1; k >= 0; -- k)
            {
                if(new_knots_m_[k] > t) new_knots_m_[k + 1] = new_knots_m_[k];
                else
                {
                    new_knots_m_[k + 1] = t;
                    // for(auto x : new_knots_m_) printf("%f ", x);
                    // printf("\n");
                    break;
                }
            }
            new_knots_length_m ++;
        }
    }

    //refine the m direction
    std::vector<std::vector<Vec3f>> new_control_points_m;
    std::vector<std::vector<float>> new_w_m;
    m = new_control_points_n[0].size();
    new_control_points_m.resize(m);
    new_w_m.resize(m);
    for(int i = 0; i < m; ++ i)
    {
        int n = new_control_points_n.size();
        new_control_points_m[i].resize(n);
        new_w_m[i].resize(n);
        for(int j = 0; j < n; ++ j)
        {
            new_control_points_m[i][j] = new_control_points_n[j][i];
            new_w_m[i][j] = new_w_n[j][i];
        }
    }
    std::vector<float> new_knots_n_;
    new_knots_n_.resize(knots_n_.size() + n_add_sum * degree_n_);
    int new_knots_length_n = 0;
    for(int i = 0; i < knots_n_.size(); ++ i) new_knots_n_[i] = knots_n_[i];
    new_knots_length_n = knots_n_.size();
    for(int ii = 0; ii < knots_n_.size() - 1; ++ ii) 
    {
        if(knots_n_[ii] == knots_n_[ii + 1]) continue;
        int ud = knots_n_addnum[ii] + 1;
        float delta = (knots_n_[ii + 1] - knots_n_[ii]) / ud;
        for(int j = 0; j < knots_n_addnum[ii]; ++ j)
        {
            float t = knots_n_[ii] + j * delta;
            int i = 0;
            for(;i + 1 < new_knots_length_n && (!(new_knots_n_[i] <= t && new_knots_n_[i + 1] > t)); ++ i);
            for(int k = 0; k < m; ++ k)
            {
                std::vector<Vec3f> tmp_controlpoints;
                auto sizet = new_control_points_m[k].size() + 1;
                tmp_controlpoints.resize(sizet);
                std::vector<float> tmp_weight;
                tmp_weight.resize(sizet);

                for(int l = 0; l <= i - degree_n_; ++ l)
                {
                    tmp_controlpoints[l] = new_control_points_m[k][l];
                    tmp_weight[l] = new_w_m[k][l];
                }
                for(int l = i - degree_n_ + 1; l <= i; ++ l)
                {
                    auto alpha = (t - new_knots_n_[l]) / (new_knots_n_[l + degree_n_] - new_knots_n_[l]);
                    tmp_controlpoints[l] = (1.0f - alpha) * new_control_points_m[k][l - 1] + alpha * new_control_points_m[k][l];
                    tmp_weight[l] = (1.0f - alpha) * new_w_m[k][l - 1] + alpha * new_w_m[k][l];
                }
                for(int l = i + 1; l < sizet; ++ l)
                {
                    tmp_controlpoints[l] = new_control_points_m[k][l - 1];
                    tmp_weight[l] = new_w_m[k][l - 1];                    
                }

                new_control_points_m[k] = tmp_controlpoints;
                new_w_m[k] = tmp_weight;
            }
            for(int k = new_knots_length_n - 1; k >= 0; -- k)
            {
                if(new_knots_n_[k] > t) new_knots_n_[k + 1] = new_knots_n_[k];
                else
                {
                    new_knots_n_[k + 1] = t;
                    break;
                }
            }
            new_knots_length_n ++;
        }
    }
    //close the m direction
    int length_n = new_knots_length_n;
    auto x_knots_n = new_knots_n_;
    for(int ii = degree_n_ + 1; ii < length_n - 1 - degree_n_; ++ ii)
    {
        for(int j = 0; j < degree_n_ - 1; ++ j)
        {
            float t = x_knots_n[ii];
            for(int k = 0; k < m; ++ k)
            {
                std::vector<Vec3f> tmp_controlpoints;
                auto sizet = new_control_points_m[k].size() + 1;
                tmp_controlpoints.resize(sizet);
                std::vector<float> tmp_weight;
                tmp_weight.resize(sizet);
                int i = 0;
                for(;i + 1 < new_knots_length_n && (!(new_knots_n_[i] <= t && new_knots_n_[i + 1] > t)); ++ i);
                for(int l = 0; l <= i - degree_n_; ++ l)
                {
                    tmp_controlpoints[l] = new_control_points_m[k][l];
                    tmp_weight[l] = new_w_m[k][l];
                }
                for(int l = i - degree_n_ + 1; l <= i; ++ l)
                {
                    auto alpha = (t - new_knots_n_[l]) / (new_knots_n_[l + degree_m_] - new_knots_n_[l]);
                    // printf("%f %f %f %f\n",new_knots_n_[l],new_knots_n_[l + degree_m_],  t, alpha);
                    tmp_controlpoints[l] = (1.0f - alpha) * new_control_points_m[k][l - 1] + alpha * new_control_points_m[k][l];
                    tmp_weight[l] = (1.0f - alpha) * new_w_m[k][l - 1] + alpha * new_w_m[k][l];
                }
                for(int l = i + 1; l < sizet; ++ l)
                {
                    tmp_controlpoints[l] = new_control_points_m[k][l - 1];
                    tmp_weight[l] = new_w_m[k][l - 1];                    
                }

                new_control_points_m[k] = tmp_controlpoints;
                new_w_m[k] = tmp_weight;
            }
            for(int k = new_knots_length_n - 1; k >= 0; -- k)
            {
                if(new_knots_n_[k] > t) new_knots_n_[k + 1] = new_knots_n_[k];
                else
                {
                    new_knots_n_[k + 1] = t;
                    break;
                }
            }
            new_knots_length_n ++;
        }
    }
    // for(int i = 0; i < new_control_points_m.size(); ++ i)
    // {
    //     for(int j = 0; j < new_control_points_m[i].size(); ++ j)
    //     {
    //         auto x = new_control_points_m[i][j];
    //         printf("%d %d (%f %f %f)\n", i, j, x.x(), x.y(), x.z());
    //     }
    // }
    //generate the mesh
    std::vector<BezierSurface> patches;
    for(int i = 0; i < new_knots_length_m - 1; ++ i)
    {
        if(new_knots_m_[i] == new_knots_m_[i + 1]) continue;
        for(int j = 0; j < new_knots_length_n - 1; ++ j)
        {
            if(new_knots_n_[j] == new_knots_n_[j + 1]) continue;
            BezierSurface patch(degree_m_ + 1, degree_n_ + 1, Vec2f(new_knots_m_[i], new_knots_m_[i + 1]), Vec2f(new_knots_n_[j], new_knots_n_[j + 1]));
            int a = 0, b = 0;
            for(int k = i - degree_m_; k <= i; ++ k)
            {
                for(int l = j - degree_n_; l <= j; ++ l)
                {
                    patch.setControlPointAndWeight(a, b, new_control_points_m[k][l] * scale + translation, new_w_m[k][l]);
                    ++ b;
                }
                ++ a;
                b = 0;
            }
            patches.push_back(patch);
        }
    }
    return std::make_shared<PatchMesh>(patches);
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
