#include <d2V_spring_particle_particle_dq2.h>

void d2V_spring_particle_particle_dq2(Eigen::Ref<Eigen::Matrix66d> H, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {
  Eigen::Vector3d delta_x = q1 - q0;
  double norm = delta_x.norm();
  double squared_norm = delta_x.squaredNorm();

  Eigen::Matrix3d term1 = (stiffness * (norm-l0) / norm) * Eigen::Matrix3d::Identity();
  Eigen::Matrix3d term2 = Eigen::Matrix3d::Constant(stiffness / squared_norm);
  Eigen::Matrix3d term3 = Eigen::Matrix3d::Constant(stiffness * (norm-l0) / (norm*squared_norm)) ;

  Eigen::Matrix3d quadrant;
  for (int i = 0; i < delta_x.size(); ++i)
    for (int j = 0; j < delta_x.size(); ++j)
      quadrant(i,j) = delta_x(i) * delta_x(j);

  H.block(0, 0, 3, 3) = H.block(3, 3, 3, 3) = quadrant.cwiseProduct(term2) + quadrant.cwiseProduct(-term3) + term1;
  H.block(3, 0, 3, 3) = H.block(0, 3, 3, 3) = quadrant.cwiseProduct(-term2) + quadrant.cwiseProduct(term3) - term1;
}