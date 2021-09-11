#include <assemble_stiffness.h>

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double k) {
  Eigen::Vector3d q0, q1;
  Eigen::Matrix66d H;
  std::vector<Eigen::Triplet<double>> coefficients;
  coefficients.reserve(q.size()*q.size());

  for (int i = 0; i < E.rows(); ++i) {
    int v0 = E(i, 0);
    int v1 = E(i, 1);
    
    q0(0) = q(v0*3);
    q0(1) = q(v0*3 + 1);
    q0(2) = q(v0*3 + 2);
    
    q1(0) = q(v1*3);
    q1(1) = q(v1*3 + 1);
    q1(2) = q(v1*3 + 2);
    
    d2V_spring_particle_particle_dq2(H, q0, q1, l0(i), k);

    int row, col;

    for (row = 0; row < 3; ++row)
      for (col = 0; col < 3; ++col)
        coefficients.emplace_back(v0 * 3 + row, v0 * 3 + col, -H(row, col));

    for (row = 0; row < 3; ++row)
      for (col = 3; col < 6; ++col)
        coefficients.emplace_back(v0 * 3 + row, v1 * 3 + col - 3, -H(row, col));

    for (row = 3; row < 6; ++row)
      for (col = 0; col < 3; ++col)
        coefficients.emplace_back(v1 * 3 + row - 3, v0 * 3 + col, -H(row, col));

    for (row = 3; row < 6; ++row)
      for (col = 3; col < 6; ++col)
        coefficients.emplace_back(v1 * 3 + row - 3, v1 * 3 + col - 3, -H(row, col));
  }
  
  K.resize(q.size(), q.size());
  K.setFromTriplets(coefficients.begin(), coefficients.end());
};
