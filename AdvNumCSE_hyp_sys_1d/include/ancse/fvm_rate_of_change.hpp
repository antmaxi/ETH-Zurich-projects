#ifndef HYPSYS1D_FVM_RATE_OF_CHANGE_HPP
#define HYPSYS1D_FVM_RATE_OF_CHANGE_HPP

#include <memory>
#include <Eigen/Dense>

#include <ancse/config.hpp>
#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/rate_of_change.hpp>
#include <ancse/simulation_time.hpp>
#include <iostream>

/// Compute the rate of change due to FVM.
/** The semidiscrete approximation of a PDE using FVM is
 *      du_i/dt = - (F_{i+0.5} - F_{i-0.5}) / dx.
 *  This computes the right hand side of the ODE.
 *
 * @tparam NumericalFlux see e.g. `CentralFlux`.
 * @tparam Reconstruction see e.g. `PWConstantReconstruction`.
 */
template <class NumericalFlux, class Reconstruction>
class FVMRateOfChange : public RateOfChange {
  public:
    FVMRateOfChange(const Grid &grid,
                    const std::shared_ptr<Model> &model,
                    const NumericalFlux &numerical_flux,
                    const Reconstruction &reconstruction)
        : grid(grid),
          model(model),
          numerical_flux(numerical_flux),
          reconstruction(reconstruction) {}

    virtual void operator()(Eigen::MatrixXd &dudt,
                            const Eigen::MatrixXd &u0) const override {
        // implement the flux loop here.

        auto n_cells = grid.n_cells;
        auto n_ghost = grid.n_ghost;
        int n_vars = model->get_nvars();
        dudt.setZero();
        double dx = grid.dx;
        //Initialize
        Eigen::VectorXd fL, fR = Eigen::VectorXd::Zero(n_vars);
        Eigen::VectorXd uL(n_vars), uR(n_vars);
        reconstruction.set(u0); // set non-changed values (buffer)
        // get flux from the last left ghost cell

        auto Ulr = reconstruction(n_ghost - 1);//|!!|..|..|
        auto Ulr1 = reconstruction(n_ghost);   //|..|!!|..|
        uL = Ulr.second;                       //|.*|..|..|
        uR = Ulr1.first;                       //|..|*.|..|
        fL = numerical_flux(uL, uR);

        for (int i = n_ghost; i < n_cells - n_ghost; ++i) {//iterating over non-ghost cells
            //                           i-1| i|i+1
            uL = Ulr1.second;          //|..|.*|..|
            Ulr1 = reconstruction(i+1);//|..|..|!!|
            uR = Ulr1.first;           //|..|..|*.|
            fR = numerical_flux(uL, uR);

            // debugging
            if (std::isnan(fR(0) + fR(1) + fR(2) + fL(0) + fL(1) + fL(2)) ) {
                std::cout << i << std::endl;
                std::cout << uL(0) << " " << uL(1) << " " << uL(2) << std::endl;
                std::cout << uR(0) << " " << uR(1) << " " << uR(2) << std::endl;
                std::cout << model->flux(uL)(0) << " " << model->flux(uL)(1) << " " << model->flux(uL)(2) << std::endl;
                std::cout << model->flux(uR)(0) << " " << model->flux(uR)(1) << " " << model->flux(uR)(2) << std::endl;
                std::cout << model->max_eigenvalue(uL) << " " << model->max_eigenvalue(uR) << std::endl;
                Eigen::VectorXd eigenval(3);
                eigenval = model->eigenvalues(uL);
                std::cout << eigenval(0) << " " << eigenval(1) << " " << eigenval(2) << std::endl;
                Eigen::VectorXd prim(3);
                prim = model->cons_to_prim(uL);
                std::cout << prim(0) << " " << prim(1) << " " << prim(2) << std::endl;
            }

            dudt.col(i) = (fL - fR) / dx;
            fL = fR;
        }
    }

  private:
    Grid grid;
    std::shared_ptr<Model> model;
    NumericalFlux numerical_flux;
    Reconstruction reconstruction;
};

std::shared_ptr<RateOfChange>
make_fvm_rate_of_change(const nlohmann::json &config,
                        const Grid &grid,
                        const std::shared_ptr<Model> &model,
                        const std::shared_ptr<SimulationTime> &simulation_time);

#endif // HYPSYS1D_FVM_RATE_OF_CHANGE_HPP
