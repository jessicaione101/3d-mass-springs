#include <fixed_point_constraints.h>
#include <algorithm>

void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {
  P.resize(q_size - 3*indices.size(), q_size);
  int row = 0;
  for (int i = 0; i < q_size/3; ++i) {
    if (std::find(indices.begin(), indices.end(), i) == indices.end()) {
      P.coeffRef(row++, i*3) = 1;
      P.coeffRef(row++, i*3 + 1) = 1;
      P.coeffRef(row++, i*3 + 2) = 1;
    }
  }
  P.makeCompressed();
}