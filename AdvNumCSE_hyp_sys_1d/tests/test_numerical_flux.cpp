#include <gtest/gtest.h>

#include <ancse/numerical_flux.hpp>

template <class NumericalFlux>
void check_consistency_euler(const NumericalFlux &nf) {

    auto model = Euler();
    int n_vars = model.get_nvars();

    Eigen::MatrixXd u(n_vars, 3);
    u.col(0) << 1, 0, 1.5;
    u.col(1) << 0.125, 0.25, 0.40;
    u.col(2) << 0.11432, -0.11432, 0.26432;

    double TOL = 1E-10;
    for (int k=0; k<u.cols(); k++) {
        ASSERT_LE( (model.flux(u.col(k)) - nf(u.col(k), u.col(k))).norm(), TOL);
    }
}

TEST(TestCentralFlux, consistency) {
    std::shared_ptr<Model> model = std::make_shared<Euler>();
    auto central_flux = CentralFlux(model);

    check_consistency_euler(central_flux);
}

TEST(TestRusanov, consistency) {
    std::shared_ptr<Model> model = std::make_shared<Euler>();
    auto rusanov = RusanovFlux(model);

    check_consistency_euler(rusanov);
}

TEST(TestLaxFriedrichs, consistency) {
    std::shared_ptr<Model> model = std::make_shared<Euler>();
    auto grid = Grid({0.0, 1.0}, 14, 2);
    auto sim_time = std::make_shared<SimulationTime>(0.0, 0.1, 0.2, 0);
    auto lf = LaxFriedrichsFlux(model, grid, sim_time);

    check_consistency_euler(lf);
}

TEST(TestRoe, consistency) {
std::shared_ptr<Model> model = std::make_shared<Euler>();
auto roe = RusanovFlux(model);

check_consistency_euler(roe);
}

TEST(TestHLLE, consistency) {
std::shared_ptr<Model> model = std::make_shared<Euler>();
auto hlle = RusanovFlux(model);

check_consistency_euler(hlle);
}

TEST(TestHLLC, consistency) {
std::shared_ptr<Model> model = std::make_shared<Euler>();
auto hllc = RusanovFlux(model);

check_consistency_euler(hllc);
}
