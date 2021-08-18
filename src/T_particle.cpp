#include <T_particle.h>
#include <iostream>

void T_particle(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, double mass) {
  T = 0.5 * mass * qdot.squaredNorm();
}
