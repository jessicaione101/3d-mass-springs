#include <assemble_forces.h>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double mass, double k) { 
  f = Eigen::VectorXd::Zero(q.size());
  Eigen::Vector3d q0, q1;
  Eigen::Vector6d forces;
  
  for (int i = 0; i < E.rows(); ++i) {
    int v0 = E(i, 0);
    int v1 = E(i, 1);
    
    q0(0) = q(v0*3);
    q0(1) = q(v0*3 + 1);
    q0(2) = q(v0*3 + 2);
    
    q1(0) = q(v1*3);
    q1(1) = q(v1*3 + 1);
    q1(2) = q(v1*3 + 2);
    
    dV_spring_particle_particle_dq(forces, q0,  q1, l0(i), k);
    
    f.segment(v0*3, 3) -= forces.segment(0, 3);
    f.segment(v1*3, 3) -= forces.segment(3, 3);
  }
};
