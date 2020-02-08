#ifndef HYPSYS1D_POLYNOMIAL_BASIS_HPP
#define HYPSYS1D_POLYNOMIAL_BASIS_HPP

#include <Eigen/Dense>
#include "grid.hpp"

/// Legendre polynomial basis
class PolynomialBasis {
  public:

    PolynomialBasis(int q, const Grid &grid) : p(q) {
        set_scaling_factor(1.0/std::sqrt(grid.dx));
    }

    PolynomialBasis(int q, double scaling_factor_) : p(q) {
        set_scaling_factor(scaling_factor_);
    }

    void set_scaling_factor(double scaling_factor_) const {
        scaling_factor = scaling_factor_;
    }

    /// Computes the Legendre polynomial basis
    /// at a given reference point xi \in [0,1]
    Eigen::VectorXd operator() (double xi) const;
    
    /// Computes the derivative of Legendre polynomial basis
    /// at a given reference point xi \in [0,1]
    Eigen::VectorXd deriv (double xi) const;

    int get_degree() const {
        return p;
    }

    double get_scaling_factor() const {
        return scaling_factor;
    }

  private:

    int p;
    mutable double scaling_factor;
};

#endif // HYPSYS1D_POLYNOMIAL_BASIS_HPP
