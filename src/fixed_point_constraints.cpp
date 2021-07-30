#include <fixed_point_constraints.h>
#include <algorithm>

void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {
  P.resize(3, q_size);
  for (int i = 0; i < indices.size(); ++i) {
    P.coeffRef(0, indices[i]) = 1;
    P.coeffRef(1, indices[i]) = 1;
    P.coeffRef(2, indices[i]) = 1;
  }
}
