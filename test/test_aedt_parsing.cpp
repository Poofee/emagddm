/**
 * @file test_aedt_parsing.cpp
 * @brief AEDT文件解析验证测试
 * @details 验证Temp.aedt文件解析的正确性
 */

#include <gtest/gtest.h>
#include "tool/maxwell_parser_impl.hpp"
#include "tool/logger_factory.hpp"
#include <filesystem>
#include <fstream>

using namespace tool;

class AedtParsingTest : public ::testing::Test {
protected:
    void SetUp() override {
        parser_ = std::make_unique<MaxwellParserImpl>();
        aedt_file_path_ = "docs/project/Temp.aedt";
        
        if (!std::filesystem::exists(aedt_file_path_)) {
            aedt_file_path_ = "../docs/project/Temp.aedt";
        }
        if (!std::filesystem::exists(aedt_file_path_)) {
            aedt_file_path_ = "../../docs/project/Temp.aedt";
        }
    }
    
    std::unique_ptr<MaxwellParserImpl> parser_;
    std::string aedt_file_path_;
};

TEST_F(AedtParsingTest, CanParseAedtFile) {
    EXPECT_TRUE(parser_->canParse(aedt_file_path_));
}

TEST_F(AedtParsingTest, ParseFileWithoutCrash) {
    ASSERT_NO_THROW({
        parser_->canParse(aedt_file_path_);
        auto materials = parser_->parseMaterials();
    }) << "Parsing Temp.aedt should not crash";
}

TEST_F(AedtParsingTest, ExtractMaterialsCount) {
    ASSERT_TRUE(parser_->canParse(aedt_file_path_));
    auto materials = parser_->parseMaterials();
    
    EXPECT_GE(materials.size(), 8) << "Should have at least 8 materials";
    
    std::cout << "Found " << materials.size() << " materials:" << std::endl;
    for (const auto& m : materials) {
        std::cout << "  - " << m["name"].get<std::string>() << std::endl;
    }
}

TEST_F(AedtParsingTest, ExtractCopperMaterial) {
    ASSERT_TRUE(parser_->canParse(aedt_file_path_));
    auto materials = parser_->parseMaterials();
    
    auto copper_it = std::find_if(materials.begin(), materials.end(),
        [](const nlohmann::json& m) {
            return m.contains("name") && m["name"].get<std::string>() == "copper";
        });
    
    ASSERT_NE(copper_it, materials.end()) << "Copper material not found";
    
    auto copper = *copper_it;
    EXPECT_NEAR(copper["conductivity"].get<double>(), 58000000.0, 1000.0);
    EXPECT_NEAR(copper["permeability"].get<double>(), 0.999991, 0.0001);
    EXPECT_NEAR(copper["mass_density"].get<double>(), 8933.0, 1.0);
}

TEST_F(AedtParsingTest, ExtractVacuumMaterial) {
    ASSERT_TRUE(parser_->canParse(aedt_file_path_));
    auto materials = parser_->parseMaterials();
    
    auto vacuum_it = std::find_if(materials.begin(), materials.end(),
        [](const nlohmann::json& m) {
            return m.contains("name") && m["name"].get<std::string>() == "vacuum";
        });
    
    ASSERT_NE(vacuum_it, materials.end()) << "Vacuum material not found";
    
    auto vacuum = *vacuum_it;
    if (vacuum.contains("permittivity")) {
        EXPECT_EQ(vacuum["permittivity"].get<double>(), 1);
    }
}

TEST_F(AedtParsingTest, ExtractNonlinearBHMaterial) {
    ASSERT_TRUE(parser_->canParse(aedt_file_path_));
    auto materials = parser_->parseMaterials();
    
    auto bh_mat_it = std::find_if(materials.begin(), materials.end(),
        [](const nlohmann::json& m) {
            return m.contains("name") && 
                   m["name"].get<std::string>().find("B27AV1400") != std::string::npos;
        });
    
    ASSERT_NE(bh_mat_it, materials.end()) << "B27AV1400 material not found";
    
    auto bh_mat = *bh_mat_it;
    
    if (bh_mat.contains("permeability_type")) {
        EXPECT_EQ(bh_mat["permeability_type"].get<std::string>(), "nonlinear");
    }
    
    if (bh_mat.contains("bh_curve")) {
        auto bh_curve = bh_mat["bh_curve"];
        EXPECT_GT(bh_curve.size(), 0) << "BH curve should have at least 1 point";
        
        if (bh_curve.size() > 0) {
            EXPECT_NEAR(bh_curve[0][0].get<double>(), 0.0, 0.001);
            EXPECT_NEAR(bh_curve[0][1].get<double>(), 0.0, 0.001);
        }
    }
}

TEST_F(AedtParsingTest, ExtractPermanentMagnetMaterial) {
    ASSERT_TRUE(parser_->canParse(aedt_file_path_));
    auto materials = parser_->parseMaterials();
    
    auto pm_it = std::find_if(materials.begin(), materials.end(),
        [](const nlohmann::json& m) {
            return m.contains("name") && 
                   m["name"].get<std::string>().find("N45UH") != std::string::npos;
        });
    
    ASSERT_NE(pm_it, materials.end()) << "N45UH permanent magnet material not found";
    
    auto pm = *pm_it;
    EXPECT_NEAR(pm["permeability"].get<double>(), 1.0739096, 0.001);
    
    ASSERT_TRUE(pm.contains("coercivity"));
    auto coercivity = pm["coercivity"];
    std::string magnitude = coercivity["magnitude"].get<std::string>();
    EXPECT_TRUE(magnitude.find("-927000") != std::string::npos);
}

TEST_F(AedtParsingTest, ExtractBoundariesCount) {
    ASSERT_TRUE(parser_->canParse(aedt_file_path_));
    parser_->parseMaterials();
    auto boundaries = parser_->parseBoundaries();
    
    EXPECT_GE(boundaries.size(), 20) << "Should have at least 20 boundaries";
    
    std::cout << "Found " << boundaries.size() << " boundaries" << std::endl;
}

TEST_F(AedtParsingTest, ExtractVectorPotentialBoundary) {
    ASSERT_TRUE(parser_->canParse(aedt_file_path_));
    parser_->parseMaterials();
    auto boundaries = parser_->parseBoundaries();
    
    auto vp_it = std::find_if(boundaries.begin(), boundaries.end(),
        [](const nlohmann::json& b) {
            return b.contains("name") && b["name"].get<std::string>() == "VectorPotential1";
        });
    
    ASSERT_NE(vp_it, boundaries.end()) << "VectorPotential1 boundary not found";
    
    auto vp = *vp_it;
    EXPECT_EQ(vp["bound_type"].get<std::string>(), "Vector Potential");
    EXPECT_EQ(vp["value"].get<std::string>(), "0");
}

TEST_F(AedtParsingTest, ExtractCoilBoundary) {
    ASSERT_TRUE(parser_->canParse(aedt_file_path_));
    parser_->parseMaterials();
    auto boundaries = parser_->parseBoundaries();
    
    auto coil_it = std::find_if(boundaries.begin(), boundaries.end(),
        [](const nlohmann::json& b) {
            return b.contains("name") && b["name"].get<std::string>() == "PhA_0";
        });
    
    ASSERT_NE(coil_it, boundaries.end()) << "PhA_0 coil boundary not found";
    
    auto coil = *coil_it;
    EXPECT_EQ(coil["bound_type"].get<std::string>(), "Coil");
    EXPECT_EQ(coil["winding"].get<int>(), 1);
    EXPECT_EQ(coil["polarity_type"].get<std::string>(), "Positive");
}

TEST_F(AedtParsingTest, ExtractExcitationsCount) {
    ASSERT_TRUE(parser_->canParse(aedt_file_path_));
    parser_->parseMaterials();
    auto excitations = parser_->parseExcitations();
    
    EXPECT_GE(excitations.size(), 10) << "Should have at least 10 coil excitations";
    
    std::cout << "Found " << excitations.size() << " excitations" << std::endl;
}

TEST_F(AedtParsingTest, ExtractSolutionSetup) {
    ASSERT_TRUE(parser_->canParse(aedt_file_path_));
    parser_->parseMaterials();
    auto setup = parser_->parseSolutionSetup();
    
    EXPECT_EQ(setup["setup_type"].get<std::string>(), "Transient");
    EXPECT_EQ(setup["stop_time"].get<std::string>(), "0.025s");
    EXPECT_EQ(setup["time_step"].get<std::string>(), "0.0001s");
    EXPECT_EQ(setup["nonlinear_residual"].get<std::string>(), "0.0001");
    EXPECT_EQ(setup["frequency"].get<std::string>(), "308.06667Hz");
}

TEST_F(AedtParsingTest, ExtractGeometry) {
    ASSERT_TRUE(parser_->canParse(aedt_file_path_));
    parser_->parseMaterials();
    auto geometry = parser_->parseGeometry();
    
    EXPECT_EQ(geometry["design_name"].get<std::string>(), "quan_ID0x45deg_20C3_temp");
    EXPECT_EQ(geometry["solution_type"].get<std::string>(), "Transient");
    EXPECT_EQ(geometry["geometry_mode"].get<std::string>(), "XY");
    EXPECT_EQ(geometry["model_depth"].get<std::string>(), "140mm");
    EXPECT_EQ(geometry["background_material"].get<std::string>(), "vacuum");
    EXPECT_EQ(geometry["units"].get<std::string>(), "mm");
}

TEST_F(AedtParsingTest, ParseAllDataIntegrity) {
    ASSERT_TRUE(parser_->canParse(aedt_file_path_));
    parser_->parseMaterials();
    auto all_data = parser_->parseAllData();
    
    EXPECT_TRUE(all_data.contains("file_info"));
    EXPECT_TRUE(all_data.contains("materials"));
    EXPECT_TRUE(all_data.contains("boundaries"));
    EXPECT_TRUE(all_data.contains("excitations"));
    EXPECT_TRUE(all_data.contains("solution_setup"));
    EXPECT_TRUE(all_data.contains("geometry"));
    
    EXPECT_FALSE(all_data["materials"].empty());
    EXPECT_FALSE(all_data["boundaries"].empty());
    
    std::cout << "\n=== Parse Summary ===" << std::endl;
    std::cout << "Materials: " << all_data["materials"].size() << std::endl;
    std::cout << "Boundaries: " << all_data["boundaries"].size() << std::endl;
    std::cout << "Excitations: " << all_data["excitations"].size() << std::endl;
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
