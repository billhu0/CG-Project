#ifndef _NURBS_H_
#define _NURBS_H_

#include <vector>
#include "core.h"
#include "geometry.h"

class NURBSCurve{
    private:
        std::vector<Vec3f> control_points_;
        std::vector<float> w_;
        std::vector<float> knots_;
        int degree_;

        int Findk(float t);
    
    public:
        NURBSCurve(std::vector<Vec3f> control_points, std::vector<float> w, std::vector<float> knots, int degree);
        NURBSCurve(int n, int degree);

        void setControlPointAndWeight(int i, Vec3f control_ponit, float w);
        void setKnot(int i, float knot);
        void setKnots(std::vector<float> knots);
        void setKnotUniform();
        std::pair<Vertex, float> evaluate(float t);
};

class NURBSSurface{
    private:
        std::vector<std::vector<Vec3f>> control_points_m_, control_points_n_;
        std::vector<std::vector<float>> w_m_, w_n_;
        std::vector<float> knots_m_, knots_n_;
        int degree_m_, degree_n_;

        // bool check_flat(Vec3f normal1, Vec3f normal2, float tolerance);
        
    public:
        NURBSSurface(int m, int n, int degreem, int degreen);

        void setControlPointAndWeight(int i, int j, Vec3f control_point, float w);
        void setKnotM(int i, float knot);
        void setKnotMUniform();
        void setKnotN(int i, float knot);
        void setKnotNUniform();
        Vertex evaluate(float u, float v);
        int findSpan(int degree, const std::vector<float> &knots, float t);
        unsigned int knotMultiplicity(const std::vector<float> &knots, float t);
        void KnotInsertU(float u, int reapeat);
        void KnotInsertV(float v, int reapeat);
        std::shared_ptr<TriangleMesh> genMesh_triangle(const Vec3f& translation, float scale);
        std::shared_ptr<PatchMesh> genMesh_patch(const Vec3f& translation, float scale);
        // Object genObject_adaptivesample(float tolerance);
};

#endif