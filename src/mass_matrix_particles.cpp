#include <mass_matrix_particles.h>

void mass_matrix_particles(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, double mass) {
  M.resize(q.size(), q.size());
  for (int i = 0; i < q.size(); ++i)
    M.coeffRef(i, i) = mass;
}
