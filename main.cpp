#include <bits/stdc++.h>
#include "CSRMatrix.hpp"

// This main implements a simple interpreter for testing and OJ:
// Input format:
// n m q
// Then q commands of the form:
//   set i j v     -> set value v at (i,j)
//   get i j       -> print value at (i,j)
//   mul k a0 .. a_{m-1} -> multiply with vector, print n values
//   slice l r     -> slice rows [l,r) and print its non-zero count and triplets
//   dense         -> print dense matrix (n lines, m values)
// Outputs are space-separated per command.

int main() {
    using sjtu::CSRMatrix;
    using sjtu::invalid_index;
    using sjtu::size_mismatch;

    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    size_t n, m, q;
    if (!(std::cin >> n >> m >> q)) return 0;
    CSRMatrix<long long> A(n, m);
    for (size_t t = 0; t < q; ++t) {
        std::string op;
        if (!(std::cin >> op)) break;
        try {
            if (op == "set") {
                size_t i, j; long long v; std::cin >> i >> j >> v;
                A.set(i, j, v);
            } else if (op == "get") {
                size_t i, j; std::cin >> i >> j;
                std::cout << A.get(i, j) << "\n";
            } else if (op == "mul") {
                // read vector of size m
                std::vector<long long> vec(m);
                for (size_t i = 0; i < m; ++i) std::cin >> vec[i];
                auto res = A * vec;
                for (size_t i = 0; i < res.size(); ++i) {
                    if (i) std::cout << ' ';
                    std::cout << res[i];
                }
                std::cout << "\n";
            } else if (op == "slice") {
                size_t l, r; std::cin >> l >> r;
                auto B = A.getRowSlice(l, r);
                std::cout << B.getRowSize() << ' ' << B.getColSize() << ' ' << B.getNonZeroCount() << "\n";
                const auto &indptr = B.getIndptr();
                const auto &indices = B.getIndices();
                const auto &data = B.getData();
                for (size_t i = 0; i < B.getRowSize(); ++i) {
                    for (size_t k = indptr[i]; k < indptr[i+1]; ++k) {
                        std::cout << i << ' ' << indices[k] << ' ' << data[k] << "\n";
                    }
                }
            } else if (op == "dense") {
                auto M = A.getMatrix();
                for (size_t i = 0; i < n; ++i) {
                    for (size_t j = 0; j < m; ++j) {
                        if (j) std::cout << ' ';
                        if (i < M.size() && j < M[i].size()) std::cout << M[i][j]; else std::cout << 0;
                    }
                    std::cout << "\n";
                }
            } else if (op == "info") {
                std::cout << A.getRowSize() << ' ' << A.getColSize() << ' ' << A.getNonZeroCount() << "\n";
            } else {
                // unknown op: consume line remainder (if any)
                std::string line; std::getline(std::cin, line);
            }
        } catch (const invalid_index &){
            std::cout << "invalid_index" << "\n";
        } catch (const size_mismatch &){
            std::cout << "size_mismatch" << "\n";
        }
    }
    return 0;
}

