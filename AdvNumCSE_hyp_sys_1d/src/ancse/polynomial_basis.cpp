#include "../../include/ancse/polynomial_basis.hpp"


/// Computes the Legendre polynomial basis
/// at a given reference point xi \in [0,1]
Eigen::VectorXd PolynomialBasis :: operator() (double xi) const
{
    Eigen::VectorXd values(p+1);
    values(0) = 1.0;
    if (p >= 1) { values(1) = xi;}
    if (p >= 2) { values(2) = 0.5*(3*xi*xi-1);}
    return scaling_factor * values;
}

/// Computes the derivative of Legendre polynomial basis
/// at a given reference point xi \in [0,1]
Eigen::VectorXd PolynomialBasis :: deriv (double xi) const
{
    Eigen::VectorXd basis_deriv(p+1);
    basis_deriv(0) = 0.0;
    if (p >= 1) { basis_deriv(1) = 1;}
    if (p >= 2) { basis_deriv(2) = 3*xi;}
    return scaling_factor * basis_deriv;
}
