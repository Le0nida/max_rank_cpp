//
// Created by leona on 02/08/2024.
//

#include <gtest/gtest.h>
#include "qnode.h"
#include "halfspace.h"
#include <vector>
#include <array>
#include <iostream>

// Helper function to create a sample Halfspace
Halfspace createHalfspace(const std::vector<double>& coeff, double known) {
    Halfspace hs;
    hs.coeff = coeff;
    hs.known = known;
    return hs;
}

// Test fixture for QNode
class QNodeTest : public ::testing::Test {
protected:
    void SetUp() override {
        std::cout << "Setting up test fixture..." << std::endl;

        // Initialize MBR for the root node
        std::vector<std::array<double, 2>> mbr = {
            std::array<double, 2>{0.0, 1.0},
            std::array<double, 2>{0.0, 1.0}
        };

        // Create root node
        root = new QNode(nullptr, mbr);
        std::cout << "Root node created." << std::endl;

        // Define the MBRs for the child nodes
        std::vector<std::array<double, 2>> mbr1 = {
            std::array<double, 2>{0.0, 0.5},
            std::array<double, 2>{0.0, 0.5}
        };
        std::vector<std::array<double, 2>> mbr2 = {
            std::array<double, 2>{0.0, 0.5},
            std::array<double, 2>{0.5, 1.0}
        };
        std::vector<std::array<double, 2>> mbr3 = {
            std::array<double, 2>{0.5, 1.0},
            std::array<double, 2>{0.0, 0.5}
        };
        std::vector<std::array<double, 2>> mbr4 = {
            std::array<double, 2>{0.5, 1.0},
            std::array<double, 2>{0.5, 1.0}
        };

        // Create child nodes and add them to the root
        auto child1 = std::make_unique<QNode>(root, mbr1);
        auto child2 = std::make_unique<QNode>(root, mbr2);
        auto child3 = std::make_unique<QNode>(root, mbr3);
        auto child4 = std::make_unique<QNode>(root, mbr4);
        child4->setNorm(false);

        root->addChild(std::move(child1));
        root->addChild(std::move(child2));
        root->addChild(std::move(child3));
        root->addChild(std::move(child4));

        // Create a few halfspaces for testing
        halfspaces = {
            createHalfspace({0.80011779, -0.21869624}, 0.5759994214752024),
            createHalfspace({0.09561655, -1.16804572}, -0.23791591035058934),
            createHalfspace({0.17754595, -0.53934412}, 0.12014104708902562),
            createHalfspace({0.61385003, 0.36678181}, 0.5323201596530454),
            createHalfspace({0.82523785, 0.32724379}, 0.6855435290914074),
            createHalfspace({0.89994708, -0.11324783}, 0.6669423488261256),
            createHalfspace({0.30242941, -0.59185904}, 0.2368170317001349),
            createHalfspace({-0.0826064, -0.77610884}, -0.09674356734052425),
            createHalfspace({-0.14718387, -0.96595473}, -0.13110948406418843),
            createHalfspace({0.9171936, 0.04861931}, 0.6684355915253815),
            createHalfspace({-0.28373698, -0.26626063}, -0.02879110680549235),
            createHalfspace({0.85476753, 0.18340153}, 0.70657481506731),
            createHalfspace({-0.14880509, -0.56388858}, -0.02355719673164358),
            createHalfspace({1.16851646, 0.63505809}, 0.7250418112802144),
            createHalfspace({0.64609332, -0.02341717}, 0.5098861313406748),
            createHalfspace({-0.46559051, -0.8240149}, -0.21091576210111407),
            createHalfspace({0.58148638, -0.44588808}, 0.4533474775747937)
        };
        std::cout << "Halfspaces created." << std::endl;

        // Create masks for testing
        masks = {
            {
                {
                    {0.0, 0.0}, {-1.0, 0.0}, {1.0, 0.0}, {0.0, -1.0}, {-1.0, -1.0}, {1.0, -1.0}, {0.0, 1.0}, {-1.0, 1.0}, {1.0, 1.0}
                },
                {
                    {1.0, 1.0, 1.0, 1.0}, {1.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 1.0}, {1.0, 0.0, 1.0, 0.0}, {1.0, 0.0, 0.0, 0.0},
                    {0.0, 0.0, 1.0, 0.0}, {0.0, 1.0, 0.0, 1.0}, {0.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 1.0}
                }
            }};
        std::cout << "Masks created." << std::endl;
    }

    void TearDown() override {
        std::cout << "Tearing down test fixture..." << std::endl;
        delete root;
        std::cout << "Root node deleted." << std::endl;
    }

    QNode* root{};
    std::vector<Halfspace> halfspaces;
    std::array<std::vector<std::vector<double>>, 2> masks;
};

// Test for isRoot
TEST_F(QNodeTest, IsRoot) {
    std::cout << "Running IsRoot test..." << std::endl;
    EXPECT_TRUE(root->isRoot());
    for (const auto& child : root->getChildren()) {
        EXPECT_FALSE(child->isRoot());
    }
    std::cout << "IsRoot test finished." << std::endl;
}

// Test for isLeaf
TEST_F(QNodeTest, IsLeaf) {
    std::cout << "Running IsLeaf test..." << std::endl;
    EXPECT_FALSE(root->isLeaf());
    for (const auto& child : root->getChildren()) {
        EXPECT_TRUE(child->isLeaf());
    }
    std::cout << "IsLeaf test finished." << std::endl;
}

// Test for getOrder
TEST_F(QNodeTest, GetOrder) {
    std::cout << "Running GetOrder test..." << std::endl;
    root->insertHalfspaces(masks, halfspaces);
    EXPECT_EQ(root->getOrder(), 0);
    for (const auto& child : root->getChildren()) {
        EXPECT_TRUE(child != nullptr); // Assicuriamoci che il child non sia nullo
        if (child->isNorm()) {
            size_t childOrder = child->getOrder();
            std::cout << "Order of child node: " << childOrder << std::endl;
            EXPECT_GT(childOrder, 0);
        }
    }
    std::cout << "GetOrder test finished." << std::endl;
}

// Test for getCovered
TEST_F(QNodeTest, GetCovered) {
    std::cout << "Running GetCovered test..." << std::endl;
    std::vector<Halfspace> covered = root->getCovered();
    std::cout << "Covered size before insertHalfspaces: " << covered.size() << std::endl;
    EXPECT_TRUE(covered.empty());

    root->insertHalfspaces(masks, halfspaces);
    covered = root->getCovered();
    std::cout << "Covered size after insertHalfspaces: " << covered.size() << std::endl;
    EXPECT_TRUE(covered.empty());
    for (const auto& child : root->getChildren()) {
        EXPECT_TRUE(child != nullptr); // Assicuriamoci che il child non sia nullo
        if (child->isNorm()) {
            covered = child->getCovered();
            std::cout << "Size of covered halfspaces: " << covered.size() << std::endl;
            EXPECT_GT(covered.size(), 0);
        }
    }
    std::cout << "GetCovered test finished." << std::endl;
}

// Test for insertHalfspaces
TEST_F(QNodeTest, InsertHalfspaces) {
    std::cout << "Running InsertHalfspaces test..." << std::endl;

    root->insertHalfspaces(masks, halfspaces);
    std::cout << "Halfspaces inserted." << std::endl;

    EXPECT_EQ(root->getChildren().size(), 4);
    for (const auto& child : root->getChildren()) {
        EXPECT_EQ(child->getChildren().size(), 0);
        if (child->isNorm() == true) {
            EXPECT_EQ(child->getCovered().size() + child->getHalfspaces().size(), halfspaces.size());
        }
    }

    std::cout << "InsertHalfspaces test finished." << std::endl;
}

// Test for setNorm and isNorm
TEST_F(QNodeTest, SetAndIsNorm) {
    std::cout << "Running SetAndIsNorm test..." << std::endl;
    root->setNorm(false);
    EXPECT_FALSE(root->isNorm());
    root->setNorm(true);
    EXPECT_TRUE(root->isNorm());
    std::cout << "SetAndIsNorm test finished." << std::endl;
}

// Test for addChild and getChildren
TEST_F(QNodeTest, AddChild) {
    std::cout << "Running AddChild test..." << std::endl;
    std::vector<std::array<double, 2>> childMbr = {
        std::array<double, 2>{0.0, 0.5},
        std::array<double, 2>{0.0, 0.5}
    };

    auto child = std::make_unique<QNode>(root, childMbr);
    root->addChild(std::move(child));

    EXPECT_FALSE(root->isLeaf());
    EXPECT_EQ(root->getChildren().size(), 5);
    std::cout << "AddChild test finished." << std::endl;
}

// Test for setHalfspaces and getHalfspaces
TEST_F(QNodeTest, SetAndGetHalfspaces) {
    std::cout << "Running SetAndGetHalfspaces test..." << std::endl;
    root->setHalfspaces(halfspaces);
    const auto& retrievedHalfspaces = root->getHalfspaces();
    std::cout << "Halfspaces size: " << retrievedHalfspaces.size() << std::endl;
    EXPECT_EQ(retrievedHalfspaces.size(), halfspaces.size());
    std::cout << "SetAndGetHalfspaces test finished." << std::endl;
}

// Test for clearHalfspaces
TEST_F(QNodeTest, ClearHalfspaces) {
    std::cout << "Running ClearHalfspaces test..." << std::endl;
    root->setHalfspaces(halfspaces);
    root->clearHalfspaces();
    EXPECT_TRUE(root->getHalfspaces().empty());
    std::cout << "ClearHalfspaces test finished." << std::endl;
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
