#include "integrator.h"
#include "utils.h"
#include <omp.h>

#include <utility>

Integrator::Integrator(std::shared_ptr<Camera> cam,
                       std::shared_ptr<Scene> scene, int spp, int max_depth)
    : camera(std::move(cam)), scene(std::move(scene)), spp(spp), max_depth(max_depth) {
}

void Integrator::render() const {
  Vec2i resolution = camera->getImage()->getResolution();
  int cnt = 0;
  Sampler sampler;
#pragma omp parallel for schedule(dynamic), shared(cnt), private(sampler)
  for (int dx = 0; dx < resolution.x(); dx++) {
#pragma omp atomic
    ++cnt;
    printf("\r%.02f%%", cnt * 100.0 / resolution.x());
    sampler.setSeed(omp_get_thread_num());
    for (int dy = 0; dy < resolution.y(); dy++) {
      Vec3f L(0, 0, 0);
      // TODO: generate #spp rays for each pixel and use Monte Carlo integration to compute radiance.
      for(int i = 0; i < spp; ++ i)
      {
        auto epsilon = sampler.get2D();
        Ray camera_ray = camera->generateRay(dx + epsilon.x(), dy + epsilon.y());
        L += radiance(camera_ray, sampler);
      }
      camera->getImage()->setPixel(dx, dy, L / float(spp));
    }
  }
}

Vec3f Integrator::radiance(Ray &ray, Sampler &sampler) const {
  Vec3f L(0, 0, 0);
  Vec3f beta(1, 1, 1);
  bool isDelta = false;
  for (int bounces = 0; ; ++bounces) {
    /// Compute radiance (direct + indirect)

    //Intersect ray with scene and store intersection
    Interaction isect;
    bool foundIntersection = scene->intersect(ray, isect);

    
    if(isect.type == Interaction::Type::LIGHT) //hit the light
    {
      if(bounces == 0) L += beta.cwiseProduct(scene->getLight()->emission(isect.pos, -ray.direction));
      // L += beta.cwiseProduct(scene->getLight()->emission(isect.pos, -ray.direction));
      break;
    }
    if(!foundIntersection || bounces >= max_depth) break;

    L += beta.cwiseProduct(directLighting(isect, sampler, isDelta));

    
    Vec3f wo = isect.wo, wi;
    float pdf = isect.material->sample(isect, sampler);
    Vec3f brdf = isect.material->evaluate(isect);
    wi = isect.wi;
    float cos_theta = isect.normal.dot(-wi);

    beta = beta.cwiseProduct(brdf) * cos_theta / pdf;
    // printf("%f %f %f\n", brdf.x() * PI, brdf.y() * PI, brdf.z() * PI);

    //RR
    // if (bounces > 3) {
    // float q = std::max((float)0.05f, 1.0f - beta.y());
    // if (sampler.get1D() < q)
    //     break;
    // beta = beta / (1.0f - q);
    // }

    ray = Ray(isect.pos, -wi);

  }
  // printf("\n");
  return L;
}

Vec3f 
Integrator::MIS(Vec3f value1, float pdf1, Vec3f value2, float pdf2) const
{
  float sum = pdf1 * pdf1 + pdf2 * pdf2;
  float w1 = pdf1 * pdf1 / sum;
  float w2 = pdf2 * pdf2 / sum;
  return value1 * w1 + value2 * w2;
}

Vec3f Integrator::directLighting(Interaction &interaction, Sampler &sampler, bool &isDelta) const {
  Vec3f L(0, 0, 0);
  // Compute direct lighting.
  float light_pdf = 0.0f;
  Vec3f light_sample_pos = scene->getLight()->sample(interaction, &light_pdf, sampler);
  Ray light_ray(interaction.pos, -interaction.wi);
  if(!scene->isShadowed(light_ray)) //sample on light
  {
    float cos_theta = interaction.normal.dot(-interaction.wi);
    float solid_angle_pdf = scene->getLight()->pdf(interaction, light_sample_pos);
    Vec3f Li = scene->getLight()->emission(light_sample_pos, -light_ray.direction);
    Vec3f fr = interaction.material->evaluate(interaction);
    L = Li.cwiseProduct(fr) * cos_theta * solid_angle_pdf / light_pdf;
  }
  if(interaction.material->isDelta())
  {
    // printf("in\n");
    float brdf_pdf = interaction.material->sample(interaction, sampler);
    Ray brdf_ray(interaction.pos, -interaction.wi);
    if(!scene->isShadowed(brdf_ray))
    {
      float cos_theta = interaction.normal.dot(-interaction.wi);
      auto Li = scene->getLight()->emission(interaction.pos, interaction.wi);
      Vec3f fr = interaction.material->evaluate(interaction);
      Vec3f L_idealspecular = Li.cwiseProduct(fr) * cos_theta / brdf_pdf;
      L = MIS(L, light_pdf, L_idealspecular, brdf_pdf);
    }
  }
  return L;
}