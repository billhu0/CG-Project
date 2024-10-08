#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

#include "camera.h"
#include "scene.h"
#include "interaction.h"

class Integrator {
   public:
    Integrator(std::shared_ptr<Camera> cam, std::shared_ptr<Scene> scene, int spp, int max_depth);
    void render() const;
    Vec3f radiance(Ray &ray, Sampler &sampler) const;

   private:
    static Vec3f MIS(const Vec3f& value1, float pdf1, const Vec3f& value2, float pdf2) ;
    Vec3f directLighting(Interaction &interaction, Sampler &sampler, bool &isDelta) const;
    std::shared_ptr<Camera> camera;
    std::shared_ptr<Scene> scene;
    int max_depth;
    int spp;
};

#endif  // INTEGRATOR_H_
