#include <ancse/dg_handler.hpp>
#include <ancse/limiters.hpp>
#include <iostream>
#include <ancse/dg_limiting.hpp>

/// build solution from DG coefficients and the basis
/// pre-evaluated at a certain point
Eigen::VectorXd DGHandler
:: build_sol(const Eigen::VectorXd& u, // n_vars*(p+1)
             const Eigen::VectorXd& basis) const
{
    Eigen::VectorXd uSol;
    uSol = Eigen::VectorXd::Zero(n_vars);
    for (int i = 0; i < n_vars; i++) {
        for (int j = 0; j < n_coeff; j++) {
            uSol(i) += u(i*n_coeff + j) * basis(j);
        }
    }
    return uSol;
}

/// build solution from DG coefficients at a given reference point
/// uses point between 0 and 1
Eigen::VectorXd DGHandler 
:: build_sol(const Eigen::VectorXd& u,
             double xi) const
{
    return build_sol(u, poly_basis(xi));
}

/// build cell average
Eigen::MatrixXd DGHandler
:: build_cell_avg (const Eigen::MatrixXd& u) const
{
    auto n_cells = u.cols();
    Eigen::MatrixXd u0 (n_vars, n_cells);
    for (int i = 0; i < n_cells; i++) {
        for (int j = 0; j < n_vars; j++) {
            u0(j, i) = u(j*n_coeff, i);
        }

    }
    return u0*poly_basis.get_scaling_factor();
}

/// build split solution uSol_m = u0 + um, uSol_p = u0 - up
/// from DG coefficients
/////////////////////////////////  ISN'T USED, EVERYTHING IN dg_limiting.cpp ////////////////////////////
std::tuple <Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>
DGHandler :: build_split_sol(const Eigen::MatrixXd& u, int n_ghost) const
{
    auto n_cells = u.cols();
    if (n_coeff > 3) {
        throw std::runtime_error(
            "Limiter not implemented for higher than 3rd order");
    }
    Eigen::MatrixXd um = Eigen::MatrixXd::Zero (n_vars, n_cells - 2*n_ghost);
    Eigen::MatrixXd up = Eigen::MatrixXd::Zero (n_vars, n_cells - 2*n_ghost);
    auto u0 = build_cell_avg(u);
    for (int i = n_ghost; i < n_cells - n_ghost; i++) {
        auto uL = build_sol(u.col(i), 0);
        auto uR = build_sol(u.col(i), 1);
        up.col(i - n_ghost) = u0.col(i) - uL;
        um.col(i - n_ghost) = uR - u0.col(i);
    }

    return {std::move(u0), std::move(um), std::move(up)}; // u0 not used
}

/// build DG coefficients from uSol_m = u0 + um, uSol_p = u0 - up
void DGHandler :: compute_limit_coeffs (Eigen::MatrixXd &u, // u to limit
                                        Eigen::MatrixXd &um, // left border
                                        Eigen::MatrixXd &up, // right border
                                        //const Limiter &limiter,
                                        int n_ghost) const
{
    std::cout << n_coeff << "!!!" << std::endl;
    if (n_coeff == 1) {
        return;
    }
    else if (n_coeff > 3) {
        throw std::runtime_error(
            "Limiter not implemented for higher than 3rd order");
    }
    else {
        Eigen::VectorXd s(n_vars);
        double scaling_factor = poly_basis.get_scaling_factor(); // for cell-averages calculation

        if (n_coeff == 2) { // um == up in this case, one unknown
            Eigen::VectorXd um_tilda(n_vars);
            Eigen::VectorXd deltaL = (u.col(n_ghost) - u.col(n_ghost - 1))*scaling_factor;
            Eigen::VectorXd deltaR;
            Eigen::VectorXd basis_coeff;
            basis_coeff = Eigen::VectorXd::Zero(n_coeff*n_vars);
            for (int i = 0; i < n_vars; i++) {
                basis_coeff(i*n_coeff + 1) = 1;
            }
            Eigen::VectorXd phi_1_m = build_sol(basis_coeff, 1); // linear part evaluated at the right border
            for (int i = n_ghost; i < u.cols() - n_ghost; i++) {
                deltaR = (u.col(i+1) - u.col(i))*scaling_factor;
                s = minmod(deltaL, deltaR);
                um_tilda = minmod(um.col(i - n_ghost), s);
                deltaL = deltaR;
                for (int coeff = 0; coeff < n_coeff - 1; coeff++) {
                    for (int var = 0; var < n_vars; var++) {
                        u(1 + coeff + var*n_coeff, i) = um_tilda(coeff*n_vars + var) /
                                                            phi_1_m(coeff*n_vars + var);
                    }
                }
                /*
                u(1, i) = um_tilda(0) / phi_1_m(0);
                u(3, i) = um_tilda(1) / phi_1_m(1);
                u(5, i) = um_tilda(2) / phi_1_m(2);
                 */
            }
    }
        else if (n_coeff == 3)
        {
                Eigen::VectorXd um_tilda(n_vars);
                Eigen::VectorXd up_tilda(n_vars);

                // preparing inverse matrix
                // linear part
                Eigen::VectorXd basis_coeff;
                basis_coeff = Eigen::VectorXd::Zero(n_coeff*n_vars);
                std::cout << n_vars << std::endl;
                for (int i = 0; i < n_vars; i++) {
                    basis_coeff(i*n_coeff + 1) = 1;
                }
                Eigen::VectorXd phi_1_m = build_sol(basis_coeff, 1); // linear part evaluated at the right border
                Eigen::VectorXd phi_1_p = build_sol(basis_coeff, 0); // linear part evaluated at the left border
                // quadratic part
                basis_coeff = Eigen::VectorXd::Zero(n_coeff*n_vars);
                for (int i = 0; i < n_vars; i++) {
                    basis_coeff(i*n_coeff + 2) = 1;
                }
                Eigen::VectorXd phi_2_m = build_sol(basis_coeff, 1); // quadratic part evaluated at the right border
                Eigen::VectorXd phi_2_p = build_sol(basis_coeff, 0); // quadratic part evaluated at the left border
                Eigen::MatrixXd Phi;
                Phi = Eigen::MatrixXd::Zero((n_coeff-1)*n_vars, (n_coeff-1)*n_vars);
                for (int i = 0; i < n_vars; i++) {
                    Phi(i, i) = phi_1_p(i);
                    Phi(i, i + n_vars) = phi_2_p(i);
                    Phi(i + n_vars, i) = phi_1_m(i);
                    Phi(i + n_vars, i + n_vars) = phi_2_m(i);
                }
                /*
                Phi << phi_1_p(0), 0, 0, phi_2_p(0), 0, 0,
                       0, phi_1_p(1), 0, 0, phi_2_p(1), 0,
                       0, 0, phi_1_p(2), 0, 0, phi_2_p(2),
                       phi_1_m(0), 0, 0, phi_2_m(0), 0, 0,
                       0, phi_1_m(1), 0, 0, phi_2_m(1), 0,
                       0, 0, phi_1_m(2), 0, 0, phi_2_m(2);
                */
                Eigen::MatrixXd Phi_inv((n_coeff-1)*n_vars, (n_coeff-1)*n_vars);
                Phi_inv = Phi.inverse();
                
                //cell average deltas
                Eigen::VectorXd deltaL = (u.col(n_ghost) - u.col(n_ghost - 1))*scaling_factor;
                Eigen::VectorXd deltaR;
                Eigen::VectorXd u_tildas((n_coeff-1)*n_vars);
                Eigen::VectorXd u_result((n_coeff-1)*n_vars);

                for (int i = n_ghost; i < u.cols() - n_ghost; i++) {
                    deltaR = (u.col(i+1) - u.col(i))*scaling_factor;
                    s = minmod(deltaL, deltaR);
                    um_tilda = minmod(um.col(i - n_ghost), s);
                    up_tilda = minmod(up.col(i - n_ghost), s);
                    deltaL = deltaR;
                    u_tildas << -up_tilda, um_tilda;
                    u_result = Phi_inv * u_tildas;

                    for (int coeff = 0; coeff < n_coeff - 1; coeff++) {
                        for (int var = 0; var < n_vars; var++) {
                            u(1 + coeff + var*n_coeff, i) = u_result(coeff*n_vars + var);
                        }
                    }
                    /*
                    u(1, i) = u_result(0); // 0,0
                    u(4, i) = u_result(1); // 0,1
                    u(7, i) = u_result(2); // 0,2
                    u(2, i) = u_result(3); // 1,0
                    u(5, i) = u_result(4); // 1,1
                    u(8, i) = u_result(5); // 1,2
                    //or
                    u(1, i) = u_result(0); // 0,0
                    u(2, i) = u_result(1); // 1,0
                     */
                }
        }
    }

}

int DGHandler :: get_n_coeff() const {
    return n_coeff;
}

int DGHandler :: get_n_vars() const {
    return n_vars;
}

double DGHandler :: get_scaling_coeff() const {
    return poly_basis.get_scaling_factor();
}