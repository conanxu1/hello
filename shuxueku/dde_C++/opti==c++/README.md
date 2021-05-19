emoji

# 🍙 功能描述

Linear DDE Solver

最优控制









# 🧂 开发环境

- g++
- eigen3
- nlopt 
- opencv
- cmake 
- make







1. sudo *apt*-get install g++ cmake make
2.  sudo *apt*-get install lib*nlopt*-dev
3. sudo apt install libopencv-dev
4. sudo apt-get install lib*eigen3*-dev 



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





# ⚽注记

条件编译

https://blog.csdn.net/weixin_42986599/article/details/106569230







