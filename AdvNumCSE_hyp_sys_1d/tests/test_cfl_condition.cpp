#include <Eigen/Dense>
#include <ancse/cfl_condition.hpp>
#include <gtest/gtest.h>

void check_cfl_condition(const CFLCondition &cfl_condition,
                         const Grid &grid,
                         double cfl_number,
                         const std::shared_ptr<Model> &model) {
    int n_ghost = grid.n_ghost;
    int n_cells = grid.n_cells;

    Eigen::MatrixXd u(3, n_cells);
    for (int i = 0; i < n_cells; ++i) {
        u(0, i) = -i * i;
        u(1, i) = -2*i * i;
        u(2, i) = -3*i * i;
    }
    double max_abs_u = (n_cells - n_ghost - 1) * (n_cells - n_ghost - 1);

    for (int i = 0; i < n_ghost; ++i) {
        for (int j = 0; j < 3; ++j) {
            u(j, i) = 100.0 * max_abs_u;
            u(j, n_cells - n_ghost + i) = 100.0 * max_abs_u;
        }
    }
    double max_v = 0.0;
    for (int i = n_ghost; i < n_cells - n_ghost; ++i) {
        max_v = std::max(max_v, model->max_eigenvalue(u.col(i)));
    }
    double dt_cfl_approx = cfl_condition(u);
    double dt_cfl_exact = cfl_number * grid.dx / max_v;

    ASSERT_DOUBLE_EQ(dt_cfl_approx, dt_cfl_exact);
}

void dg_check_cfl_condition(
                            const CFLCondition &cfl_condition,
                            const Grid &grid,
                         double cfl_number,
                            const DGHandler &dg_handler,
                         const std::shared_ptr<Model> &model) {
    int n_ghost = grid.n_ghost;
    int n_cells = grid.n_cells;
    int n_coeff = dg_handler.get_n_coeff();
    Eigen::MatrixXd u(3*n_coeff, n_cells);
    for (int i = 0; i < n_cells; ++i) {
        for (int j = 0; j < 3*n_coeff; j++) {
            if (j % n_coeff == 0) {
                u(j,i) = -(j % n_coeff + 1) * i * i;
            } else {
                u(j,i) = 100 * i * i;
            }
        }

    }
    double max_abs_u = (n_cells - n_ghost - 1) * (n_cells - n_ghost - 1);

    for (int i = 0; i < n_ghost; ++i) {
        for (int j = 0; j < 3*n_coeff; ++j) {
            u(j, i) = 100.0 * max_abs_u;
            u(j, n_cells - n_ghost + i) = 100.0 * max_abs_u;
        }
    }
    auto u0 = dg_handler.build_cell_avg(u);
    double max_v = 0.0;
    for (int i = n_ghost; i < n_cells - n_ghost; ++i) {
        max_v = std::max(max_v, model->max_eigenvalue(u0.col(i)));
    }
    double dt_cfl_approx = cfl_condition(u);
    double dt_cfl_exact = cfl_number * grid.dx / max_v;

    ASSERT_DOUBLE_EQ(dt_cfl_approx, dt_cfl_exact);
}

TEST(CFLCondition, Example) {
auto n_ghost = 2;
auto n_cells = 10 + 2 * n_ghost;
auto grid = Grid({0.9, 1.0}, n_cells, n_ghost);
std::shared_ptr<Model> model = std::make_shared<Euler>();
double cfl_number = 0.2;

auto cfl_condition = StandardCFLCondition(grid, model, cfl_number);
check_cfl_condition(cfl_condition, grid, cfl_number, model);
}

TEST(DG_CFLCondition, Example) {
auto n_ghost = 2;
auto n_cells = 10 + 2 * n_ghost;
auto grid = Grid({0.9, 1.0}, n_cells, n_ghost);
std::shared_ptr<Model> model = std::make_shared<Euler>();
double cfl_number = 0.2;

auto dg_handler = DGHandler(model, PolynomialBasis(2, grid));
auto cfl_condition = DG_CFLCondition(grid, model, dg_handler, cfl_number);
dg_check_cfl_condition(cfl_condition, grid, cfl_number, dg_handler, model);
}