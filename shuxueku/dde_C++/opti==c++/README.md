emoji

# 🍙 功能描述

Linear DDE Solver

最优控制









# 🧂 开发环境

- g++
- eigen3
- nlopt



# 🍬 **项目结构**



# 🍼 **DEMO**

采用cmake编译build_by_cmake
opti==c++$ ./*sh "-DUSE_EIGEN=1 -DSMALL_SCALE=1"
使用EIGEN库，求解小规模问题不优化



opti==c++$ ./*sh "-DUSE_EIGEN=1 -DLARGE_SCALE=1"
使用EIGEN库，求解大规模问题。采用数组保存数据，计算过程中复制数据给MatrixXd，或者可以考虑使用openblas库进行处理。










# 🍺 **作者列表**

conanxu1







# **🍪历史版本**













