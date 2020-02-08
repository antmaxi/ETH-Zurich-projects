#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <iostream>
#include <ancse/model.hpp>

TEST(TestModel, Euler) {
    std::shared_ptr<Model> model = std::make_shared<Euler>();
    Eigen::VectorXd u(3);
    u << 2, 3, 4; // rho, moment, E

    auto euler = std::dynamic_pointer_cast<Euler>(model);
    double gamma = euler->get_gamma();
    Eigen::VectorXd u_prim(3);
    u_prim << 2, 1.5, (gamma - 1)*(4 - 0.5*3*3/2.0); // rho, v, p
    Eigen::VectorXd u_to_check(3);
    u_to_check = euler->cons_to_prim(u);
    // check cons_to_prim
    for (int i = 0; i < 3; i++) {
        ASSERT_DOUBLE_EQ(u_to_check(i), u_prim(i));
    }
    // check prim_to_cons
    for (int i = 0; i < 3; i++) {
        ASSERT_DOUBLE_EQ(euler->prim_to_cons(u_to_check)(i), u(i));
    }
    // check flux
    Eigen::VectorXd flux_to_check(3);
    flux_to_check << 3, 3*3/2.0 + u_prim(2), (4 + u_prim(2))*3/2.0;
    for (int i = 0; i < 3; i++) {
        ASSERT_DOUBLE_EQ(euler->flux(u)(i), flux_to_check(i));
    }
    // check getCH
    Eigen::VectorXd CH_to_check(2);
    CH_to_check << std::sqrt(gamma*u_prim(2)/u(0)), (u(2) + u_prim(2))/u(0);
    for (int i = 0; i < 2; i++) {
        ASSERT_DOUBLE_EQ(euler->getCH(u)(i), CH_to_check(i));
    }
    // check eigenvalues
    Eigen::VectorXd eigenvalues(3);
    double c = CH_to_check(0);
    double v = u_prim(1);
    eigenvalues << v - c, v, v + c;
    for (int i = 0; i < 3; i++) {
        ASSERT_DOUBLE_EQ(euler->eigenvalues(u)(i), eigenvalues(i));
    }
    // check eigenvalues
    Eigen::MatrixXd eigenvectors(3, 3);
    double h = CH_to_check(1);
    eigenvectors << 1, v-c, h-v*c,
                    1,   v, v*v/2.0,
                    1, v+c, h+v*c;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ASSERT_DOUBLE_EQ(euler->eigenvectors(u)(i, j), eigenvectors(i, j));
        }
    }
    //check max eigenvalue
    double max_eigenvalue = std::max(std::max(std::abs(eigenvalues(0)), std::abs(eigenvalues(1))), std::abs(eigenvalues(2)));
    ASSERT_DOUBLE_EQ(euler->max_eigenvalue(u), max_eigenvalue);
}

