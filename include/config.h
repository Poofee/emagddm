// 项目配置头文件
// 定义编译期开关宏，用于控制MPI和OpenMP功能

#ifndef FE_EM_CONFIG_H
#define FE_EM_CONFIG_H

// MPI编译开关
// 0=关闭MPI，1=开启MPI
#ifndef USE_MPI
#define USE_MPI 0
#endif

// OpenMP编译开关
// 0=关闭OpenMP，1=开启OpenMP
#ifndef USE_OPENMP
#define USE_OPENMP 0
#endif

// MPI根目录（由CMake自动检测并定义）
#ifndef MPI_ROOT
#define MPI_ROOT ""
#endif

// 编译期模式检测宏
#if USE_MPI && USE_OPENMP
#define HYBRID_MODE_MPI_OMP 1
#define SERIAL_MODE 0
#define PURE_OMP_MODE 0
#define PURE_MPI_MODE 0
#elif USE_MPI
#define HYBRID_MODE_MPI_OMP 0
#define SERIAL_MODE 0
#define PURE_OMP_MODE 0
#define PURE_MPI_MODE 1
#elif USE_OPENMP
#define HYBRID_MODE_MPI_OMP 0
#define SERIAL_MODE 0
#define PURE_OMP_MODE 1
#define PURE_MPI_MODE 0
#else
#define HYBRID_MODE_MPI_OMP 0
#define SERIAL_MODE 1
#define PURE_OMP_MODE 0
#define PURE_MPI_MODE 0
#endif

// 模式描述字符串
#if HYBRID_MODE_MPI_OMP
#define COMPILATION_MODE "MPI+OpenMP混合模式"
#elif PURE_MPI_MODE
#define COMPILATION_MODE "纯MPI模式"
#elif PURE_OMP_MODE
#define COMPILATION_MODE "纯OpenMP模式"
#else
#define COMPILATION_MODE "串行模式"
#endif

#endif // FE_EM_CONFIG_H