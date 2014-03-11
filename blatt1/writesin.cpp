#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

int main()
{
    double x = .0;

    std::cout<<"# Plotting sin fuction in [0,2pi] with delta=0.01pi.";
    std::cout<<std::scientific;

    while (x < 2* M_PI)
    {
        std::cout<<x<<"  "<<std::sin(x)<<'\n';
        x += .01* M_PI;
    };

    return 0;
}

