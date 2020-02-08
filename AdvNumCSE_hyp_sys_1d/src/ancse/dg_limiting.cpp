#include <ancse/dg_limiting.hpp>

#include <iostream>
#include <ancse/config.hpp>
#include <ancse/boundary_condition.hpp>


/// applies the DG limiting procedure to the solution coefficients
/// for all cells.
template <class Limiter>
void DGLimiting <Limiter> :: operator()(Eigen::MatrixXd &u) const
{
    Eigen::MatrixXd u0, um, up;
    std::tie(u0, um, up) = dg_handler.build_split_sol(u, n_ghost);
    //dg_handler.compute_limit_coeffs(u, um, up, //limiter,
    //        n_ghost); // get limited coefficients

    int n_vars = dg_handler.get_n_vars();
    int n_coeff = dg_handler.get_n_coeff();
    Eigen::VectorXd s(n_vars);
    Eigen::VectorXd um_tilda(n_vars);
    Eigen::VectorXd up_tilda(n_vars);
    Eigen::VectorXd basis_coeff;
    basis_coeff = Eigen::VectorXd::Zero(n_coeff*n_vars);
    // basis for linear part:
    for (int i = 0; i < n_vars; i++) {
        basis_coeff(i*n_coeff + 1) = 1;
    }

    if (n_coeff == 2) { // um == up in this case, one unknown
        Eigen::VectorXd phi_1_m = dg_handler.build_sol(basis_coeff, 1); // linear part evaluated at the right border
        for (int i = n_ghost; i < u.cols() - n_ghost; i++) {
            std::tie(um_tilda, up_tilda) = (*this)(u0.col(i - 1), u0.col(i), u0.col(i+1),
                                                         um.col(i - n_ghost), up.col(i - n_ghost));
            for (int coeff = 0; coeff < n_coeff - 1; coeff++) {
                for (int var = 0; var < n_vars; var++) {
                    u(1 + coeff + var * n_coeff, i) = um_tilda(coeff * n_vars + var) /
                                                      phi_1_m(coeff * n_vars + var);
                }
            }
        }
    } else if (n_coeff == 3)
    {
        ////////////////////// preparing inverse matrix
        // first linear part:
        Eigen::VectorXd phi_1_m = dg_handler.build_sol(basis_coeff, 1); // linear part evaluated at the right border
        Eigen::VectorXd phi_1_p = dg_handler.build_sol(basis_coeff, 0); // linear part evaluated at the left border
        // quadratic part
        basis_coeff = Eigen::VectorXd::Zero(n_coeff*n_vars);
        for (int i = 0; i < n_vars; i++) {
            basis_coeff(i*n_coeff + 2) = 1;
        }
        Eigen::VectorXd phi_2_m = dg_handler.build_sol(basis_coeff, 1); // quadratic part evaluated at the right border
        Eigen::VectorXd phi_2_p = dg_handler.build_sol(basis_coeff, 0); // quadratic part evaluated at the left border
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

        Eigen::VectorXd u_tildas((n_coeff-1)*n_vars);
        Eigen::VectorXd u_result((n_coeff-1)*n_vars);

        for (int i = n_ghost; i < u.cols() - n_ghost; i++) {
            std::tie(um_tilda, up_tilda) = (*this)(u0.col(i - 1), u0.col(i), u0.col(i+1),
                                                   um.col(i - n_ghost), up.col(i - n_ghost));
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
            /*
            u(1, i) = um_tilda(0) / phi_1_m(0);
            u(3, i) = um_tilda(1) / phi_1_m(1);
            u(5, i) = um_tilda(2) / phi_1_m(2);

        }
    }
    */
}

/// DG limiting procedure
/** uc0 : cell average of cell 'i' (inscaled)
 *  um0 : cell average of cell 'i-1'
 *  up0 : cell average of cell 'i+1'
 *  uSol_{@right_face_of_cell_i} = uc0 + um
 *  uSol_{@left_face_of_cell_i}  = uc0 - up
*/
template <class Limiter>
inline std::pair<Eigen::VectorXd, Eigen::VectorXd>
DGLimiting <Limiter> :: operator()(const Eigen::VectorXd &um0,
                                   const Eigen::VectorXd &uc0,
                                   const Eigen::VectorXd &up0,
                                   const Eigen::VectorXd &um,
                                   const Eigen::VectorXd &up) const
{
    double scaling_factor = dg_handler.get_scaling_coeff();
    Eigen::VectorXd deltaL = (uc0 - um0) * scaling_factor;
    Eigen::VectorXd deltaR = (up0 - uc0) * scaling_factor;
    Eigen::VectorXd um_tilda(uc0.size());
    Eigen::VectorXd up_tilda(uc0.size());
    for (int i = 0; i < dg_handler.get_n_vars(); i++) {
        um_tilda(i) = limiter(deltaL(i), deltaR(i), um(i));
        up_tilda(i) = limiter(deltaL(i), deltaR(i), up(i));
    }

    return {std::move(um_tilda), std::move(up_tilda)};
}


#define REGISTER_DG_LIMITER(token, LimiterType, limiter)      \
    if (config["dg_limiter"] == token) {                      \
        return std::make_shared<DGLimiting<LimiterType>> (    \
            grid, dg_handler, limiter);   \
    }

std::shared_ptr<Limiting>
make_dg_limiting(const nlohmann::json &config,
                 const Grid &grid,
                 const DGHandler &dg_handler)
{
    REGISTER_DG_LIMITER("vanleer", VanLeer, VanLeer{})

    REGISTER_DG_LIMITER("shu", Shu, Shu(grid.dx))

    throw std::runtime_error(fmt::format(
        "Unknown DG limiter. [{}]", std::string(config["dg_limiter"])));
}

#undef REGISTER_DG_LIMITER
