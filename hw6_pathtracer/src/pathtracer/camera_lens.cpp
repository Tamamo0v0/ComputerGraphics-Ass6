#include "camera.h"

#include <iostream>
#include <sstream>
#include <fstream>

#include "CGL/misc.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"

using std::cout;
using std::endl;
using std::max;
using std::min;
using std::ifstream;
using std::ofstream;

namespace CGL {

using Collada::CameraInfo;

Ray Camera::generate_ray_for_thin_lens(double x, double y, double rndR, double rndTheta) const {

  // TODO Assignment 7: Part 4
  // compute position and direction of ray from the input sensor sample coordinate.
  // Note: use rndR and rndTheta to uniformly sample a unit disk.
  
  double w_bias = tan(hFov * M_PI / 180 / 2);
  double h_bias = tan(vFov * M_PI / 180 / 2);

  double x_pix = (x - 0.5) * 2 * w_bias;
  double y_pix = (y - 0.5) * 2 * h_bias;


  Vector3D dir = {x_pix,y_pix,-1};
  Vector3D pFocus = focalDistance * dir;


  Vector3D pLens = {lensRadius * sqrt(rndR) * cos(rndTheta), lensRadius * sqrt(rndR) * sin(rndTheta), 0};

  Vector3D ray_pos = c2w * pLens + pos;
  Vector3D ray_dir = c2w * (pFocus - pLens).unit();

  ray_dir.normalize();
  
  Ray ray(ray_pos, ray_dir);
  ray.min_t = nClip; 
  ray.max_t = fClip;

  return ray;

}


} // namespace CGL
