#!/bin/bash
# MPI/OpenMP编译开关验证脚本
# 测试四种编译组合的编译和运行情况

set -e

echo "=== MPI/OpenMP编译开关验证 ==="
echo ""

# 清理之前的构建
if [ -d "build" ]; then
    echo "清理之前的构建..."
    rm -rf build
fi

# 测试组合1: USE_MPI=0, USE_OPENMP=0 (串行模式)
echo "=== 测试组合1: 串行模式 (USE_MPI=0, USE_OPENMP=0) ==="
mkdir -p build
cd build
cmake .. -DUSE_MPI=OFF -DUSE_OPENMP=OFF
make -j4
./test/test_omp_basic
./test/test_mpi_basic
cd ..
rm -rf build
echo "✅ 组合1测试通过"
echo ""

# 测试组合2: USE_MPI=0, USE_OPENMP=1 (纯OpenMP模式)
echo "=== 测试组合2: 纯OpenMP模式 (USE_MPI=0, USE_OPENMP=1) ==="
mkdir -p build
cd build
cmake .. -DUSE_MPI=OFF -DUSE_OPENMP=ON
make -j4
./test/test_omp_basic
cd ..
rm -rf build
echo "✅ 组合2测试通过"
echo ""

# 测试组合3: USE_MPI=1, USE_OPENMP=0 (纯MPI模式)
echo "=== 测试组合3: 纯MPI模式 (USE_MPI=1, USE_OPENMP=0) ==="
mkdir -p build
cd build
cmake .. -DUSE_MPI=ON -DUSE_OPENMP=OFF
make -j4
# MPI测试需要特殊运行方式
mpirun -np 2 ./test/test_mpi_basic
cd ..
rm -rf build
echo "✅ 组合3测试通过"
echo ""

# 测试组合4: USE_MPI=1, USE_OPENMP=1 (混合模式)
echo "=== 测试组合4: 混合模式 (USE_MPI=1, USE_OPENMP=1) ==="
mkdir -p build
cd build
cmake .. -DUSE_MPI=ON -DUSE_OPENMP=ON
make -j4
# 混合模式测试
mpirun -np 2 ./test/test_mpi_basic
./test/test_omp_basic
cd ..
rm -rf build
echo "✅ 组合4测试通过"
echo ""

echo "🎉 所有编译组合测试通过！"
echo "MPI/OpenMP编译开关和无感适配功能验证完成"