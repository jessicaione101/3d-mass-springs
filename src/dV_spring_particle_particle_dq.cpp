#include <dV_spring_particle_particle_dq.h>

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {
  Eigen::Vector3d delta_x = q1 - q0;
  double norm = delta_x.norm();
  f << -delta_x, delta_x;
  f = f * stiffness * (norm - l0) / norm;
}
