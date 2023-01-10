#include "bsdf.h"
#include "utils.h"

#include <utility>

const float default_delta = 1e-7;

IdealDiffusion::IdealDiffusion(const Vec3f &color) : color(color) { }

Vec3f IdealDiffusion::evaluate(Interaction &interaction) const {
    return color * INV_PI;
}

float IdealDiffusion::pdf(Interaction &interaction) const {
    float cos_theta = interaction.normal.dot(-interaction.wi);
    return cos_theta * INV_PI;
}

float IdealDiffusion::sample(Interaction &interaction, Sampler &sampler) const {
    Vec2f epsilon = sampler.get2D();
    float r = sqrtf(epsilon.x()), phi = 2.0f * PI * epsilon.y();
    float x = r * cosf(phi), y = r * sinf(phi);
    float z = sqrtf(std::max(0.0f, 1.0f - x * x - y * y));
    Vec3f local_wi(-x, -y, -z);

    float rotate_theta = acosf(Vec3f(0, 0, 1).dot(interaction.normal));
    Vec3f rotate_axis = Vec3f(0, 0, 1).cross(interaction.normal).normalized();
    Eigen::AngleAxisf rotationVec(rotate_theta, rotate_axis);
    Eigen::Matrix3f rotation_matrix = rotationVec.toRotationMatrix();
    interaction.wi = (rotation_matrix * local_wi).normalized();
    float cos_theta = interaction.normal.dot(-interaction.wi);
    return cos_theta * INV_PI;
}

// return whether the bsdf is perfect transparent or perfect reflection
bool IdealDiffusion::isDelta() const {
    return false;
}

// Ideal Specular
IdealSpecular::IdealSpecular(const Vec3f &color) : color(color) { }

Vec3f IdealSpecular::evaluate(Interaction &interaction) const {
    return color * pdf(interaction);
}

float IdealSpecular::pdf(Interaction &interaction) const {
    Vec3f reflect = (-interaction.wo + 2.0f * interaction.wo.dot(interaction.normal) * interaction.normal).normalized();
    if ((reflect + interaction.wi).norm() < default_delta) {
        return 1.0f;
    } else {
        return 0.0f;
    }
}

float IdealSpecular::sample(Interaction &interaction, Sampler &sampler) const {
    interaction.wi = (interaction.wo - 2.0f * interaction.wo.dot(interaction.normal) * interaction.normal).normalized();
    return 1.0f;
}

bool IdealSpecular::isDelta() const {
    return true;
}