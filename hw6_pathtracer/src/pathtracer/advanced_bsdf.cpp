#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>

#include "application/visual_debugger.h"

using std::max;
using std::min;
using std::swap;

namespace CGL {

// Mirror BSDF //

Vector3D MirrorBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D MirrorBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

  // TODO Assignment 7: Part 1
  // Implement MirrorBSDF
  return Vector3D();
}

void MirrorBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Mirror BSDF"))
  {
    DragDouble3("Reflectance", &reflectance[0], 0.005);
    ImGui::TreePop();
  }
}

// Microfacet BSDF //

double MicrofacetBSDF::G(const Vector3D wo, const Vector3D wi) {
  return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D h) {
  // TODO Assignment 7: Part 2
  // Compute Beckmann normal distribution function (NDF) here.
  // You will need the roughness alpha.

  double theta = getTheta(h.unit());
  double tan_theta = tan(theta);
  double cos_theta = cos(theta);

  double res = exp(-((tan_theta * tan_theta) / (alpha * alpha))) 
                / (M_PI * alpha * alpha * cos_theta * cos_theta * cos_theta * cos_theta);


  return res;
}

Vector3D MicrofacetBSDF::F(const Vector3D wi) {
  // TODO Assignment 7: Part 2
  // Compute Fresnel term for reflection on dielectric-conductor interface.
  // You will need both eta and etaK, both of which are Vector3D.

  double theta = getTheta(wi.unit());
  double cos_theta = cos(theta);

  Vector3D Rs = ((eta * eta + k * k) - 2 * eta * cos_theta + cos_theta * cos_theta)
                /((eta * eta + k * k) + 2 * eta * cos_theta + cos_theta * cos_theta);
  Vector3D Rp = ((eta * eta + k * k) * cos_theta * cos_theta - 2 * eta * cos_theta + 1)
                /((eta * eta + k * k) * cos_theta * cos_theta + 2 * eta * cos_theta + 1);

  Vector3D res = (Rs + Rp) / 2;

  return res;
}

Vector3D MicrofacetBSDF::f(const Vector3D wo, const Vector3D wi) {
  // TODO Assignment 7: Part 2
  // Implement microfacet model here.

  Vector3D n = {0,0,1.};

  if(dot(wo,n) < 0 || dot(wi,n) < 0){
    return Vector3D(0,0,0);
  }

  Vector3D h = (wo + wi) / 2;
  h.normalize();

  Vector3D res = (F(wi) * G(wo,wi) * D(h)) / (4 * dot(n,wo) * dot(n,wi));

  return res;
}

Vector3D MicrofacetBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // TODO Assignment 7: Part 2
  // *Importance* sample Beckmann normal distribution function (NDF) here.
  // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
  //       and return the sampled BRDF value.

  //default
  //*wi = cosineHemisphereSampler.get_sample(pdf);
  //return MicrofacetBSDF::f(wo, *wi);

  //importance
  double r1 = random_uniform();
  double r2 = random_uniform();

  double theta = atan(sqrt(-(alpha * alpha) * log(1 - r1))); 
  double phi = 2 * M_PI * r2;

  Vector3D h = {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
  h.normalize();
  
  *wi = dot(wo,h) * h - wo + dot(wo,h) * h;
  wi->normalize();

  Vector3D n = {0,0,1.};

  if(dot(*wi,n) < 0){
    *pdf = 0;
    return Vector3D(0,0,0);
  }

  double p_theta = (2 * sin(theta) / (alpha * alpha * cos(theta) * cos(theta) * cos(theta)))
                    * exp(-((tan(theta) * tan(theta)) / (alpha * alpha)));
  double p_phi = 1.0 / (M_PI * 2);

  double p_h = p_theta * p_phi / sin(theta);
  double p_wi = p_h / (4 * dot(*wi,h));

  *pdf = p_wi;

  return MicrofacetBSDF::f(wo, *wi);
}

void MicrofacetBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Micofacet BSDF"))
  {
    DragDouble3("eta", &eta[0], 0.005);
    DragDouble3("K", &k[0], 0.005);
    DragDouble("alpha", &alpha, 0.005);
    ImGui::TreePop();
  }
}

// Refraction BSDF //

Vector3D RefractionBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D RefractionBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // TODO Assignment 7: Part 1
  // Implement RefractionBSDF
  return Vector3D();
}

void RefractionBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Refraction BSDF"))
  {
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("ior", &ior, 0.005);
    ImGui::TreePop();
  }
}

// Glass BSDF //

Vector3D GlassBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D GlassBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

  // TODO Assignment 7: Part 1
  // Compute Fresnel coefficient and either reflect or refract based on it.

  // compute Fresnel coefficient and use it as the probability of reflection
  // - Fundamentals of Computer Graphics page 305
  return Vector3D();
}

void GlassBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Refraction BSDF"))
  {
    DragDouble3("Reflectance", &reflectance[0], 0.005);
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("ior", &ior, 0.005);
    ImGui::TreePop();
  }
}

void BSDF::reflect(const Vector3D wo, Vector3D* wi) {

  // TODO Assignment 7: Part 1
  // Implement reflection of wo about normal (0,0,1) and store result in wi.


}

bool BSDF::refract(const Vector3D wo, Vector3D* wi, double ior) {

  // TODO Assignment 7: Part 1
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.

  return true;

}

} // namespace CGL
