#include<pro.h>

void test()
{
Mat coef = (Mat_<float>(3,1) << 1,-2,1);
Mat roots;
solvePoly(coef, roots);
cout << "Roots: channels = " << roots.channels() << " , values = " << roots << ".";

}
