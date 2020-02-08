#include <ancse/model.hpp>

#include <iostream>
#include <fmt/format.h>


///------------------///
/// Euler equations  ///
///------------------///
// u - in conservative by default
Eigen::VectorXd Euler::flux(const Eigen::VectorXd &u) const
{
    Eigen::VectorXd flux(3);
    Eigen::VectorXd u_prim = cons_to_prim(u); // rho, v, p
    flux << u(1),
            u(1)*u_prim(1) + u_prim(2),
            (u(2) + u_prim(2))*u_prim(1);
    return flux;
}

// get velocity c and enthalpy H
Eigen::VectorXd Euler::getCH(const Eigen::VectorXd &u) const
{
    Eigen::VectorXd u_prim = cons_to_prim(u);
    Eigen::VectorXd ch(2);
    ch << std::sqrt(gamma*u_prim(2)/u(0)), // c
          (u(2)+u_prim(2))/u(0);   // H
    return ch;
}

Eigen::VectorXd Euler::eigenvalues(const Eigen::VectorXd &u) const
{
    Eigen::VectorXd eigenval(3);
    double c = getCH(u)(0);
    double v = u(1)/u(0);
    eigenval << v-c, v, v+c;
    //std::cout << v << " -- " << c << std::endl;
    return eigenval;
}

Eigen::MatrixXd Euler::eigenvectors(const Eigen::VectorXd &u) const
{
    Eigen::MatrixXd eigenvec(3, 3);
    Eigen::VectorXd ch(2);
    double v = u(1)/u(0);
    ch = getCH(u);
    eigenvec << 1, v - ch(0), ch(1) - v*ch(0),
                1,         v,           v*v/2,
                1, v + ch(0), ch(1) + v*ch(0);
    return eigenvec;
}

double Euler::max_eigenvalue(const Eigen::VectorXd &u) const
{
    Eigen::VectorXd eigenval(3);
    eigenval = eigenvalues(u);
    return std::max(std::abs(eigenval(0)), std::abs(eigenval(2)));
}

// get (rho, v, p) from (rho, momentum, energy)
Eigen::VectorXd Euler::cons_to_prim(const Eigen::VectorXd &u_cons) const
{
    Eigen::VectorXd u_prim(3);
    u_prim << std::max(u_cons(0), 1e-10),
              u_cons(1) / u_cons(0),//std::max(0.0,
               std::max(1e-10, (gamma-1)*(u_cons(2) - 0.5*u_cons(1)*u_cons(1)/u_cons(0)));
    return u_prim;
}

// get (rho, momentum, energy) from (rho, v, p)
Eigen::VectorXd Euler::prim_to_cons(const Eigen::VectorXd &u_prim) const
{
    Eigen::VectorXd u_cons(3);
    u_cons << u_prim(0),
              u_prim(1)*u_prim(0),
              u_prim(2)/(gamma - 1) + 0.5*u_prim(0)*u_prim(1)*u_prim(1);
    return u_cons;
}


#define REGISTER_MODEL(token, ModelType)      \
    if (config["model"] == (token)) {         \
        return std::make_shared<ModelType>(); \
    }

std::shared_ptr<Model> make_model (const nlohmann::json &config)
{
    REGISTER_MODEL("burgers", Burgers)
    REGISTER_MODEL("euler", Euler)
    // implement and register your models here

    throw std::runtime_error(
        fmt::format("Unknown model. {}", std::string(config["flux"])));
}
