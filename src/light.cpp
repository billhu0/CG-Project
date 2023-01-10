#include "light.h"

#include <utility>
#include "utils.h"

constexpr uint32_t SAMPLE_NUM = 16;

Light::Light(Vec3f pos, Vec3f color) : position(std::move(pos)), radiance(std::move(color)) { }

SquareAreaLight::SquareAreaLight(const Vec3f &pos, const Vec3f &color, const Vec2f &size)
    : Light(pos, color), size(size), normal(Vec3f(0, -1, 0)) {
    Vec3f v1, v2, v3, v4;
    v1 = pos + Vec3f(size.x() / 2, 0.f, -size.y() / 2);
    v2 = pos + Vec3f(-size.x() / 2, 0.f, -size.y() / 2);
    v3 = pos + Vec3f(-size.x() / 2, 0.f, size.y() / 2);
    v4 = pos + Vec3f(size.x() / 2, 0.f, size.y() / 2);
    light_mesh = TriangleMesh({v1, v2, v3, v4}, {Vec3f(0, -1, 0)}, {0, 1, 2, 0, 2, 3}, {0, 0, 0, 0, 0, 0});
}

Vec3f SquareAreaLight::emission(const Vec3f &pos, const Vec3f &dir) const {
    if (dir.dot(normal) > 0.0f) {
        return radiance;
    } else {
        return Vec3f::Zero();
    }
}

float SquareAreaLight::pdf(const Interaction &interaction, Vec3f pos) {
    if (abs((pos - position).dot(Vec3f(1, 0, 0))) > size.x() / 2.0f ||
        abs((pos - position).dot(Vec3f(0, 0, 1))) > size.y() / 2.0f)
        return 0;
    float cos_theta1 = abs(normal.dot(interaction.wi));
    float dis = (pos - interaction.pos).norm();
    return cos_theta1 / (dis * dis);
}

Vec3f SquareAreaLight::sample(Interaction &interaction, float *pdf, Sampler &sampler) const {
    auto epsilon = sampler.get2D();
    Vec3f tangentx(1, 0, 0), tangentz(0, 0, 1);
    auto pos = position + (epsilon.x() - 0.5f) * size.x() * tangentx + (epsilon.y() - 0.5f) * size.y() * tangentz;
    interaction.wi = (interaction.pos - pos).normalized();
    *pdf = 1.0f / (size.x() * size.y());
    return pos;
}

bool SquareAreaLight::intersect(Ray &ray, Interaction &interaction) const {
    if (light_mesh.intersect(ray, interaction)) {
        interaction.type = Interaction::Type::LIGHT;
        return true;
    }
    return false;
}
