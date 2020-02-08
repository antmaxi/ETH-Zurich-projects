#ifndef HYPSYS1D_NUMERICAL_FLUX_HPP
#define HYPSYS1D_NUMERICAL_FLUX_HPP

#include <memory>
#include <iostream>

#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/simulation_time.hpp>

/// Central flux.
/** This flux works does not depend on the model.
 * It is also unconditionally a bad choice.
 */
class CentralFlux {
  public:
    // Note: the interface for creating fluxes will give you access
    //       to the following:
    //         - model
    //         - grid
    //         - shared_ptr to simulation_time
    //       Therefore, try to only use a subset of those three in your
    //       constructors.
    explicit CentralFlux(const std::shared_ptr<Model> &model)
        : model(model) {}

    /// Compute the numerical flux given the left and right trace.
    Eigen::VectorXd operator()(const Eigen::VectorXd &uL,
                               const Eigen::VectorXd &uR) const
    {
        auto fL = model->flux(uL);
        auto fR = model->flux(uR);

        return 0.5 * (fL + fR);
    }

  private:
    std::shared_ptr<Model> model;
};

class RusanovFlux {
public:
    explicit RusanovFlux(const std::shared_ptr<Model> &model)
            : model(model) {}

    /// Compute the numerical flux given the left and right trace.
    Eigen::VectorXd operator()(const Eigen::VectorXd &uL,
                               const Eigen::VectorXd &uR) const
    {
        auto fL = model->flux(uL);
        auto fR = model->flux(uR);

        return 0.5 * (fL + fR) - 0.5 * std::max(model->max_eigenvalue(uL), model->max_eigenvalue(uR))*(uR - uL);
    }

private:
    std::shared_ptr<Model> model;
};

class LaxFriedrichsFlux {
public:
    explicit LaxFriedrichsFlux(const std::shared_ptr<Model> &model,
                               const Grid &grid,
                               std::shared_ptr<SimulationTime> simulation_time)
            : model(model),
              grid(grid),
              simulation_time(std::move(simulation_time)) {}

    /// Compute the numerical flux given the left and right trace.
    Eigen::VectorXd operator()(const Eigen::VectorXd &uL,
                               const Eigen::VectorXd &uR) const
    {
        double dx = grid.dx;
        double dt = simulation_time->dt;
        auto fL = model->flux(uL);
        auto fR = model->flux(uR);

        return 0.5 * ((fL + fR) - dx / dt * (uR - uL));
    }

private:
    std::shared_ptr<Model> model;
    std::shared_ptr<SimulationTime> simulation_time;
    Grid grid;
};

class RoeFlux {
        public:
        explicit RoeFlux(const std::shared_ptr<Model> &model)
        : model(model) {}

        /// Compute the numerical flux given the left and right trace.
        Eigen::VectorXd operator()(const Eigen::VectorXd &uL,
        const Eigen::VectorXd &uR) const
        {
            auto euler = std::dynamic_pointer_cast<Euler>(model);
            double gamma = euler->get_gamma();

            auto fL = model->flux(uL);
            auto fR = model->flux(uR);

            double rho_sqrL = std::sqrt(uL(0));
            double rho_sqrR = std::sqrt(uR(0));
            Eigen::VectorXd uL_pr = model->cons_to_prim(uL);
            Eigen::VectorXd uR_pr = model->cons_to_prim(uR);

            // get Roe's variables
            double v_hat = (rho_sqrL * uL_pr(1) + rho_sqrR * uR_pr(1)) /
                           (rho_sqrL + rho_sqrR);
            double v_hat2 = v_hat * v_hat;
            double H_hat = (rho_sqrL * (uL(2) + uL_pr(2)) / uL(0)
                            + rho_sqrR * (uR(2) + uR_pr(2)) / uR(0)) /
                           (rho_sqrL + rho_sqrR);
            /*
            Eigen::Matrix2Xd RoeMatrix(3, 3);
            RoeMatrix << 0, 1, 0,
                    0.5*(gamma-3)*v_hat2, (3-gamma)*v_hat, gamma-1,
                    (0.5*(gamma-1)*v_hat2-H_hat)*v_hat, H_hat - (gamma-1)*v_hat2, gamma*v_hat;
            */
            double c_hat = std::sqrt((gamma-1)*(H_hat - 0.5*v_hat2));
            double c_hat2 = c_hat * c_hat;
            // Eigenvalues of Roe's matrix
            Eigen::MatrixXd Lambda(3,3);
            Lambda << std::abs(v_hat - c_hat), 0,         0,
                      0,     std::abs(v_hat),         0,
                      0,         0, std::abs(v_hat + c_hat);
            Eigen::MatrixXd R(3,3);
            // Eigenvectors of Roe's matrix
            R <<             1,     1,             1,
                 v_hat - c_hat, v_hat, v_hat + c_hat,
                 H_hat - v_hat*c_hat, 0.5*v_hat2, H_hat + v_hat*c_hat;
            Eigen::MatrixXd R_inv(3,3);
            R_inv << 0.25*(gamma-1)*v_hat2/c_hat2 + 0.5*v_hat/c_hat,
                    -0.5*(gamma-1)*v_hat/c_hat2 - 0.5/c_hat,
                    0.5*(gamma-1)/c_hat2,
                    
                    1 - 0.5*(gamma-1)*v_hat2/c_hat2, 
                    (gamma-1)*v_hat/c_hat2,
                    (gamma-1)/c_hat2,
                    
                    0.25*(gamma-1)*v_hat2/c_hat2 - 0.5*v_hat/c_hat,
                    -0.5*(gamma-1)*v_hat/c_hat2 + 0.5/c_hat,
                    0.5*(gamma-1)/c_hat2;                    
            return 0.5 * (fL + fR) - 0.5 *R*Lambda*R_inv*(uR - uL);
        }

        private:
        std::shared_ptr<Model> model;
};

class HLLEFlux {
public:
    explicit HLLEFlux(const std::shared_ptr<Model> &model)
            : model(model) {}

    /// Compute the numerical flux given the left and right trace.
    Eigen::VectorXd operator()(const Eigen::VectorXd &uL,
                               const Eigen::VectorXd &uR) const
    {

        auto euler = std::dynamic_pointer_cast<Euler>(model);
        double gamma = euler->get_gamma();

        double rho_sqrL = std::sqrt(uL(0));
        double rho_sqrR = std::sqrt(uR(0));
        Eigen::VectorXd uL_pr = model->cons_to_prim(uL);
        Eigen::VectorXd uR_pr = model->cons_to_prim(uR);

        // get Roe's variables
        double v_hat = (rho_sqrL * uL_pr(1) + rho_sqrR * uR_pr(1)) /
                       (rho_sqrL + rho_sqrR);
        double v_hat2 = v_hat * v_hat;
        double H_hat = (rho_sqrL * (uL(2) + uL_pr(2)) / uL(0)
                        + rho_sqrR * (uR(2) + uR_pr(2)) / uR(0)) /
                       (rho_sqrL + rho_sqrR);


        double c_hat = std::sqrt((gamma-1)*(H_hat - 0.5*v_hat2));

        double cL = std::sqrt(gamma * uL_pr(2) / uL_pr(0));
        double cR = std::sqrt(gamma * uR_pr(2) / uR_pr(0));

        // HLL wave speeds
        double sL = std::min(std::min(std::min(v_hat - c_hat, uL_pr(1) - cL),
                                      std::min(v_hat, uL_pr(1))),
                                      std::min(v_hat + c_hat, uL_pr(1) + cL));
        double sR = std::max(std::max(std::max(v_hat - c_hat, uR_pr(1) - cR),
                                      std::max(v_hat, uR_pr(1))),
                                      std::max(v_hat + c_hat, uR_pr(1) + cR));
        if (sL >= 0) {
            return model->flux(uL);
        } else if (sR <= 0) {
            return model->flux(uR);
        } else {
            auto fL = model->flux(uL);
            auto fR = model->flux(uR);
            return (sR*fL - sL*fR + sL*sR*(uR-uL))/(sR-sL);
        }
}

private:
    std::shared_ptr<Model> model;
};

class HLLCFlux {
public:
    explicit HLLCFlux(const std::shared_ptr<Model> &model)
            : model(model) {}

    /// Compute the numerical flux given the left and right trace.
    Eigen::VectorXd operator()(const Eigen::VectorXd &uL,
                               const Eigen::VectorXd &uR) const {

        double rho_sqrL = std::sqrt(uL(0));
        double rho_sqrR = std::sqrt(uR(0));
        Eigen::VectorXd uL_pr = model->cons_to_prim(uL);
        Eigen::VectorXd uR_pr = model->cons_to_prim(uR);

        // get Roe's variables
        double v_hat = (rho_sqrL * uL_pr(1) + rho_sqrR * uR_pr(1)) /
                       (rho_sqrL + rho_sqrR);
        double v_hat2 = v_hat * v_hat;
        double H_hat = (rho_sqrL * (uL(2) + uL_pr(2)) / uL(0)
                        + rho_sqrR * (uR(2) + uR_pr(2)) / uR(0)) /
                       (rho_sqrL + rho_sqrR);
        auto euler = std::dynamic_pointer_cast<Euler>(model);
        double gamma = euler->get_gamma();

        double c_hat = std::sqrt((gamma - 1) * (H_hat - 0.5 * v_hat2));

        double cL = std::sqrt(gamma * uL_pr(2) / uL_pr(0));
        double cR = std::sqrt(gamma * uR_pr(2) / uR_pr(0));

        // HLL wave speeds
        double sL = std::min(std::min(std::min(v_hat - c_hat, uL_pr(1) - cL),
                                      std::min(v_hat, uL_pr(1))),
                             std::min(v_hat + c_hat, uL_pr(1) + cL));
        double sR = std::max(std::max(std::max(v_hat - c_hat, uR_pr(1) - cR),
                                      std::max(v_hat, uR_pr(1))),
                             std::max(v_hat + c_hat, uR_pr(1) + cR));
        if (sL >= 0) {
            return model->flux(uL);
        } else if (sR <= 0) {
            return model->flux(uR);
        } else {
            // middle speed
            double sM = (uR_pr(0) * uR_pr(1) * (sR - uR_pr(1))
                         - uL_pr(0) * uL_pr(1) * (sL - uL_pr(1))
                         - (uR_pr(2) - uL_pr(2))) /
                        (uR_pr(0) * (sR - uR_pr(1)) - uL_pr(0) * (sL - uL_pr(1)));
            double p_star = uR_pr(2) + uR_pr(0) * (uR_pr(1) - sM) * (uR_pr(1) - sR);
            if (sM > 0) {
                double rhoL = uL_pr(0) * (uL_pr(1) - sL) / (sM - sL);
                Eigen::VectorXd uM(3);
                uM << rhoL, rhoL*sM, p_star/(gamma-1) + 0.5*rhoL*sM*sM;
                return model->flux(uL) + sL*(uM-uL);
            } else {
                double rhoR = uR_pr(0) * (uR_pr(1) - sR) / (sM - sR);
                Eigen::VectorXd uM(3);
                uM << rhoR, rhoR*sM, p_star/(gamma-1) + 0.5*rhoR*sM*sM;
                return model->flux(uR) + sR*(uM-uR);
            }
        }
    }

private:
    std::shared_ptr<Model> model;
};

#endif // HYPSYS1D_NUMERICAL_FLUX_HPP
