#include <fixed_point_constraints.h>
#include <algorithm>

void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {
  P.resize(3*(q_size - indices.size()), q_size);
  P.setZero();
  for (int i = 0; i < indices.size(); ++i) {
    int col = P.cols() - indices[i]*3 - 1;
    P.coeffRef(i*3,     col)   = 1;
    P.coeffRef(i*3 + 1, col-1) = 1;
    P.coeffRef(i*3 + 2, col-2) = 1;
  }
}
