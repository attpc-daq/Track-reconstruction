# Track-reconstruction
2D径迹重构-Hough Transform

2D直线重构Hough变换具体可参考以下网站进行理解：

https://blog.csdn.net/nichongben/article/details/79029966

https://blog.csdn.net/weixin_42192493/article/details/105765597


其中main.cpp是用Hough变换重构径迹的一个示例。由于直线重构精度与划分bin的大小有关，但随着bin划分越小，Hough变换计算量越大，耗时越长。因此在main.cpp中采用了先粗略bin下的Hough变换重建径迹，将所得的参数$\left (\theta,r  \right )$ 作为函数TF1预拟合参数，拟合得到精确值。
