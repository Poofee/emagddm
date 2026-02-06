#include <iostream>
#include <nlohmann/json.hpp>

int main() {
    // 测试nlohmann/json库是否可用
    nlohmann::json j;
    j["test"] = "Hello, World!";
    j["number"] = 42;
    
    std::cout << "JSON库测试成功: " << j.dump(4) << std::endl;
    
    return 0;
}