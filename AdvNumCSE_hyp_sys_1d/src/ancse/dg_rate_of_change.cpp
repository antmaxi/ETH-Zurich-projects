#include <ancse/dg_rate_of_change.hpp>

#include <Eigen/Dense>
#include <ancse/config.hpp>
#include <ancse/polynomial_basis.hpp>
#include <ancse/dg_handler.hpp>
#include <ancse/numerical_flux.hpp>
#include <fmt/format.h>
#include <iostream>

/// DG numerical flux term
template <class NumericalFlux>
void DGRateOfChange<NumericalFlux>
:: eval_numerical_flux (Eigen::MatrixXd &dudt, //n_vars*(p+1) X n_cells
                        const Eigen::MatrixXd &u0) const //n_vars*(p+1) X n_cells
{
    // implement the loop for DG numerical flux term.

    int n_vars = model->get_nvars();
    int n_cells = u0.cols();
    int n_coeff = poly_basis.get_degree() + 1;
    Eigen::VectorXd poly_basisL(n_vars);
    poly_basisL = poly_basis(0); // basis evaluated on the left side
    Eigen::VectorXd poly_basisR(n_vars);
    poly_basisR = poly_basis(1); //                 on the right side

    for (int j = grid.n_ghost; j <= n_cells - grid.n_ghost; j++) {// iterating over interfaces
        Eigen::VectorXd uL(n_vars); // U at the left side of interface
        Eigen::VectorXd uR(n_vars); //   at the right
        uL = dg_handler.build_sol(u0.col(j - 1), poly_basisR);
        uR = dg_handler.build_sol(u0.col(j), poly_basisL);
        if (0) {
            std::cout << "|" << j << "|" << std::endl;
            for (int a = 0; a < n_vars * n_coeff; a++) {
                std::cout << u0.col(j - 1)(a) << " ";
            }
            std::cout << std::endl;
            for (int a = 0; a < n_vars * n_coeff; a++) {
                std::cout << u0.col(j)(a) << " ";
            }
            std::cout << std::endl;
            std::cout << uL(0) << " " << uL(1) << " " << uL(2) << " !! " << uL.size() << std::endl;
            std::cout << uR(0) << " " << uR(1) << " " << uR(2) << " !! " << uR.size() << std::endl;
            std::cout << n_vars << std::endl;
        }
        auto flux = numerical_flux(uL, uR);
        for (int i = 0; i < n_vars; i++) {
            for (int l = 0; l < n_coeff; l++) {
                dudt(i * n_coeff + l, j) += flux(i) * poly_basisL(l);
                dudt(i * n_coeff + l, j - 1) -= flux(i) * poly_basisR(l);
            }
        }

    }

}

/// DG volume integral term
template <class NumericalFlux>
void DGRateOfChange<NumericalFlux>
:: eval_volume_integral(Eigen::MatrixXd &dudt, //n_vars*(p+1) X n_cells
                        const Eigen::MatrixXd &u0) const //n_vars*(p+1) X n_cells
{
    // implement the loop for DG volume integral.

    int n_vars = model->get_nvars();
    int n_cells = u0.cols();
    int n_coeff = poly_basis.get_degree() + 1;
    int n_ghost = grid.n_ghost;
    Eigen::VectorXd u(n_vars);
    if (n_coeff >= 2) {
        for (int k = 0; k < quad_weights.size(); k++) {// over points in quadrature
            auto basis = poly_basis(quad_points(k));
            auto deriv_basis = poly_basis.deriv(quad_points(k));
            for (int i = n_ghost; i < n_cells - n_ghost; i++) {// iterating over cells
                u = dg_handler.build_sol(u0.col(i), basis);
                for (int j = 0; j < n_vars; j++) {
                    auto flux = model->flux(u);
                    for (int l = 0; l < n_coeff; l++) {
                        dudt(j * n_coeff + l, i) += quad_weights(k) * flux(j) * deriv_basis(l);
                    }
                }
            }
        }
    }

}
#define REGISTER_NUMERICAL_FLUX(token, FluxType, flux)          \
    if (config["flux"] == (token)) {                            \
        return std::make_shared< DGRateOfChange<FluxType> >(    \
            grid, model, flux, poly_basis, dg_handler);                     \
    }


std::shared_ptr<RateOfChange> make_dg_rate_of_change(
    const nlohmann::json &config,
    const Grid &grid,
    const std::shared_ptr<Model> &model,
    const PolynomialBasis &poly_basis,
    const DGHandler &dg_handler,
    const std::shared_ptr<SimulationTime> &simulation_time)
{
    // Register the other numerical fluxes.
    REGISTER_NUMERICAL_FLUX("central_flux", CentralFlux, CentralFlux(model))

    // Register the other numerical fluxes.
    REGISTER_NUMERICAL_FLUX("rusanov", RusanovFlux, RusanovFlux((model)));
    REGISTER_NUMERICAL_FLUX("lax_friedrichs", LaxFriedrichsFlux,
                            LaxFriedrichsFlux(model, grid, simulation_time));
    REGISTER_NUMERICAL_FLUX("roe", RoeFlux, RoeFlux((model)));
    REGISTER_NUMERICAL_FLUX("hlle", HLLEFlux, HLLEFlux((model)));
    REGISTER_NUMERICAL_FLUX("hllc", HLLCFlux, HLLCFlux(model));
    throw std::runtime_error(
        fmt::format("Unknown numerical flux. {}",
                    std::string(config["flux"])));
}

#undef REGISTER_NUMERICAL_FLUX
