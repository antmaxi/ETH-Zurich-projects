#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <ancse/reconstruction.hpp>


TEST(TestPWConstant, Example) {
    auto rc = PWConstantReconstruction{};

    Eigen::VectorXd ua(3), ub(3);
    ua << 1.0, 0, 1.50;
    ub << 0.1, 0, 0.15;


    auto [uL, uR] = rc(ua, ub);

    ASSERT_EQ(uL, ua);
    ASSERT_EQ(uR, ub);
}

TEST(TestMinMod, Example) {
ASSERT_DOUBLE_EQ(minmod(-1.0, -2.0), -1.0);
ASSERT_DOUBLE_EQ(minmod(2.0, 1.0), 1.0);

ASSERT_DOUBLE_EQ(minmod(-1.0, 3.0), 0.0);

ASSERT_DOUBLE_EQ(minmod(0.0, -1.0), 0.0);
ASSERT_DOUBLE_EQ(minmod(1.0, 0.0), 0.0);
}


TEST(TestConsPWLinear, Example) {

auto rc = ConsPWLinearReconstruction{MinMod{}};

Eigen::VectorXd ua(3);
ua << 1.5, 1.5, 1.5;
Eigen::VectorXd ub(3);
ub << 2.0, 2.0, 2.0;
Eigen::VectorXd uc(3);
uc << 3.0, 3.0, 3.0;
auto [uL, uR] = rc(ua, ub, uc);

for  (int i = 0; i < 3; i++) {
    ASSERT_DOUBLE_EQ(uL(i), 1.75);
    ASSERT_DOUBLE_EQ(uR(i), 2.25);
}

}

TEST(TestPrimPWLinear, Example) {
std::shared_ptr<Model> model = std::make_shared<Euler>();
auto rc = PrimPWLinearReconstruction{MinMod{}, model};

Eigen::VectorXd ua(3);
ua << 1.5, 1.5, 1.5;
Eigen::VectorXd ub(3);
ub << 2.0, 2.0, 2.0;
Eigen::VectorXd uc(3);
uc << 3.0, 3.0, 3.0;

auto [uL, uR] = rc(ua, ub, uc);

for  (int i = 0; i < 3; i++) {
    ASSERT_DOUBLE_EQ(uL(i), 1.75);
    ASSERT_DOUBLE_EQ(uR(i), 2.25);
}
}