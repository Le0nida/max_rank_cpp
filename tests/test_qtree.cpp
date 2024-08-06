//
// Created by leona on 01/08/2024.
//

#include <gtest/gtest.h>
#include "halfspace.h"
#include "qtree.h"
#include "qnode.h"
#include <vector>
#include <iostream>

// Helper function to create a sample Halfspace
HalfSpace createHalfspace(const std::vector<double>& coeff, double known)
{
    Point dummyPoint({0.0, 0.0}); // Assuming a dummy point for construction
    Eigen::VectorXd coeff_eigen = Eigen::Map<const Eigen::VectorXd>(coeff.data(), coeff.size());
    return {dummyPoint, coeff_eigen, known};
}

// Test fixture for QTree
class QTreeTest : public ::testing::Test {
protected:
    void SetUp() override {
        std::cout << "Setting up test fixture..." << std::endl;

        // Initialize QTree
        dims = 2;
        maxhsnode = 4;
        qtree = std::make_unique<QTree>(dims, maxhsnode);

        // Create a few halfspaces for testing
        halfspaces = {
            createHalfspace({0.8, -0.2}, 0.57),
            createHalfspace({0.1, -1.2}, -0.24),
            createHalfspace({0.2, -0.5}, 0.12),
            createHalfspace({0.6, 0.3}, 0.53),
            createHalfspace({0.8, 0.3}, 0.68),
            createHalfspace({0.9, -0.1}, 0.67),
            createHalfspace({0.3, -0.6}, 0.24),
            createHalfspace({-0.1, -0.8}, -0.10),
            createHalfspace({-0.1, -1.0}, -0.13),
            createHalfspace({0.9, 0.0}, 0.67)
        };

        std::cout << "Halfspaces created." << std::endl;
    }

    void TearDown() override {
        std::cout << "Tearing down test fixture..." << std::endl;
    }

    int dims{};
    int maxhsnode{};
    std::unique_ptr<QTree> qtree;
    std::vector<HalfSpace> halfspaces;
};

// Test for inserting halfspaces
TEST_F(QTreeTest, InsertHalfspace) {
    std::cout << "Running InsertHalfspace test..." << std::endl;
    qtree->inserthalfspaces(halfspaces);
    auto leaves = qtree->getleaves();
    EXPECT_GT(leaves.size(), 0);
    std::cout << "InsertHalfspace test finished." << std::endl;
}

// Test for getting leaves
TEST_F(QTreeTest, GetLeaves) {
    std::cout << "Running GetLeaves test..." << std::endl;
    auto leaves = qtree->getleaves();
    EXPECT_GT(leaves.size(), 0);
    std::cout << "GetLeaves test finished." << std::endl;
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
