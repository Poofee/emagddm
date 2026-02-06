/**
 * @file test_tool_modules.cpp
 * @brief 基础工具层模块单元测试程序
 * @details 测试operation_logger、version_manager等功能
 * @author AI Developer
 * @date 2026-02-06
 * @version 1.0
 */

#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <memory>
#include <cstdio>
#include <algorithm>
#include <chrono>
#include <cstdlib>

#include "tool/operation_logger.hpp"
#include "tool/version_manager.hpp"
#include "tool/file_utils.hpp"

class OperationLoggerTest : public ::testing::Test {
protected:
    void SetUp() override {
    }
    
    void TearDown() override {
    }
};

TEST_F(OperationLoggerTest, LogCreate) {
    tool::OperationLogger logger;
    
    uint64_t log_id = logger.log_create("Material", "Copper", "Created copper material");
    
    EXPECT_NE(log_id, 0);
    
    auto logs = logger.get_all_logs();
    EXPECT_EQ(logs.size(), 1);
    EXPECT_EQ(logs[0].target_type, "Material");
    EXPECT_EQ(logs[0].target_id, "Copper");
}

TEST_F(OperationLoggerTest, LogRead) {
    tool::OperationLogger logger;
    
    uint64_t log_id = logger.log_read("Geometry", "MotorCore", "Read geometry data");
    
    EXPECT_NE(log_id, 0);
    
    auto logs = logger.get_logs_by_type(tool::OperationType::READ);
    EXPECT_EQ(logs.size(), 1);
}

TEST_F(OperationLoggerTest, LogWrite) {
    tool::OperationLogger logger;
    
    uint64_t log_id = logger.log_write("Boundary", "InnerBoundary", "Updated boundary condition");
    
    EXPECT_NE(log_id, 0);
}

TEST_F(OperationLoggerTest, LogModify) {
    tool::OperationLogger logger;
    
    uint64_t log_id = logger.log_modify(
        "Material", 
        "Copper", 
        "Updated conductivity",
        "5.8e7",
        "5.9e7"
    );
    
    EXPECT_NE(log_id, 0);
    
    auto logs = logger.get_logs_by_type(tool::OperationType::MODIFY);
    ASSERT_EQ(logs.size(), 1);
    EXPECT_EQ(logs[0].before_value, "5.8e7");
    EXPECT_EQ(logs[0].after_value, "5.9e7");
}

TEST_F(OperationLoggerTest, LogDelete) {
    tool::OperationLogger logger;
    
    uint64_t log_id = logger.log_delete("Excitation", "Coil1", "Removed coil excitation");
    
    EXPECT_NE(log_id, 0);
    
    auto logs = logger.get_logs_by_type(tool::OperationType::DELETE);
    EXPECT_EQ(logs.size(), 1);
}

TEST_F(OperationLoggerTest, LogImport) {
    tool::OperationLogger logger;
    
    uint64_t log_id = logger.log_import("motor.aedt", "Geometry", "Imported geometry from Maxwell");
    
    EXPECT_NE(log_id, 0);
    
    auto logs = logger.get_logs_by_type(tool::OperationType::IMPORT);
    EXPECT_EQ(logs.size(), 1);
}

TEST_F(OperationLoggerTest, LogExport) {
    tool::OperationLogger logger;
    
    uint64_t log_id = logger.log_export(
        "Mesh", 
        "MotorMesh", 
        "exported_mesh.vtk", 
        "Exported mesh to VTK format"
    );
    
    EXPECT_NE(log_id, 0);
    
    auto logs = logger.get_logs_by_type(tool::OperationType::EXPORT);
    EXPECT_EQ(logs.size(), 1);
}

TEST_F(OperationLoggerTest, GetLogsByTarget) {
    tool::OperationLogger logger;
    
    logger.log_create("Material", "Copper", "Created");
    logger.log_modify("Material", "Copper", "Modified");
    logger.log_delete("Material", "Aluminum", "Deleted");
    
    auto logs = logger.get_logs_by_target("Material", "Copper");
    EXPECT_EQ(logs.size(), 2);
}

TEST_F(OperationLoggerTest, UpdateStatus) {
    tool::OperationLogger logger;
    
    uint64_t log_id = logger.log_create("Material", "Test", "Test");
    
    bool result = logger.update_operation_status(log_id, tool::OperationStatus::FAILED, "Test error");
    EXPECT_TRUE(result);
    
    auto log = logger.get_log(log_id);
    ASSERT_TRUE(log.has_value());
    EXPECT_EQ(log->status, tool::OperationStatus::FAILED);
    EXPECT_TRUE(log->additional_info.find("Test error") != std::string::npos);
}

TEST_F(OperationLoggerTest, ExportToFile) {
    tool::OperationLogger logger;
    
    logger.log_create("Material", "Copper", "Created copper");
    
    std::string export_file = "test_operation_log.xml";
    bool result = logger.export_to_file(export_file);
    EXPECT_TRUE(result);
    
    std::ifstream file(export_file);
    ASSERT_TRUE(file.is_open());
    
    std::string content((std::istreambuf_iterator<char>(file)),
                        std::istreambuf_iterator<char>());
    
    EXPECT_TRUE(content.find("OperationLogs") != std::string::npos);
    EXPECT_TRUE(content.find("Copper") != std::string::npos);
    
    file.close();
    std::remove(export_file.c_str());
}

TEST_F(OperationLoggerTest, OperatorName) {
    tool::OperationLogger logger;
    logger.set_operator_name("TestOperator");
    
    logger.log_create("Material", "TestMat", "Created");
    
    auto logs = logger.get_all_logs();
    ASSERT_EQ(logs.size(), 1);
    EXPECT_EQ(logs[0].operator_name, "TestOperator");
}

class VersionDiffTest : public ::testing::Test {
protected:
    void SetUp() override {
    }
    
    void TearDown() override {
    }
};

TEST_F(VersionDiffTest, AddEntry) {
    tool::VersionDiff::DiffEntry entry;
    entry.data_type = "Material";
    entry.entity_id = "Copper";
    entry.field_name = "conductivity";
    entry.old_value = "5.8e7";
    entry.new_value = "5.9e7";
    
    tool::VersionDiff diff;
    diff.add_entry(entry);
    
    EXPECT_EQ(diff.get_entry_count(), 1);
    EXPECT_EQ(diff.get_modified_count(), 1);
}

TEST_F(VersionDiffTest, AddNewEntity) {
    tool::VersionDiff diff;
    diff.add_new_entity("Material", "NewMaterial");
    
    EXPECT_EQ(diff.get_entry_count(), 1);
    EXPECT_EQ(diff.get_new_count(), 1);
}

TEST_F(VersionDiffTest, AddDeletedEntity) {
    tool::VersionDiff diff;
    diff.add_deleted_entity("Boundary", "OldBoundary");
    
    EXPECT_EQ(diff.get_entry_count(), 1);
    EXPECT_EQ(diff.get_deleted_count(), 1);
}

TEST_F(VersionDiffTest, ToString) {
    tool::VersionDiff diff;
    diff.add_new_entity("Material", "Mat1");
    diff.add_deleted_entity("Boundary", "Bnd1");
    diff.add_modified_field("Material", "Mat2", "value", "old", "new");
    
    std::string str = diff.to_string();
    EXPECT_TRUE(str.find("VersionDiff") != std::string::npos);
    EXPECT_TRUE(str.find("3") != std::string::npos);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
