#include <map>
#include <valarray>

#include <iostream>
#include <fstream>

#define _USE_MATH_DEFINES
#include <cmath>

#include "coplib.hpp"

// initialize with default paramters
static std::map<std::string, double> params {
        {"stepsize", 0.01},
        {"resolution", 1},
        {"t_init", 0.0}, // initial conditions
        {"t_end", 50},
        {"y1_init", 0},
        {"y2_init", 0},
        {"y3_init", 4},
        {"y4_init", 2},
        {"omega", 1}
};

static double omega;

std::valarray<double> dglfunc_doppelpendel(double /* t */, const std::valarray<double>& y)
{
    std::valarray<double> ret(4);

    double sin_dtheta = std::sin(y[0] - y[1]), cos_dtheta = std::cos(y[0] - y[1]);

    double A =  (y[2]*y[3]*sin_dtheta) / (1 + sin_dtheta*sin_dtheta);
    double B =  (y[2]*y[2] + 2*y[3]*y[3] - 2*y[2]*y[3]*cos_dtheta) / 
                 (1 + sin_dtheta*sin_dtheta)*(1 + sin_dtheta*sin_dtheta) *
                 sin_dtheta*cos_dtheta;

    // dottheta_1 
    ret[0] = (y[2] - y[3]*cos_dtheta) / (1 + sin_dtheta*sin_dtheta);
    
    // dottheta_2
    ret[1] = (2*y[3] -y[2]*cos_dtheta) / (1 + sin_dtheta*sin_dtheta);
    
    // dot-p-tilde_1
    ret[2] = -A + B - 2*omega*omega*std::sin(y[0]);
    
    // dot-p-tilde_2
    ret[3] =  A - B - omega*omega*std::sin(y[1]);


    return ret;
}

double energy(std::valarray<double> y)
{
    return (1 - std::cos(y[0])) + y[1]*y[1]/2;
}

int main()
{
    std::ifstream parmfile("03_doublependulum.prm");
    if (!parmfile)
    {
        perror("Failed to open parameter file");
        std::cerr<<"Using default parameters\n";
    }
    else
    {
        cop::readParams(params, parmfile);
        std::clog<<"Parameters loaded.\n";
    }




	// initialize solver
    double timestep = params["stepsize"], t_end = params["t_end"];
    omega = params["omega"];
    std::valarray<double> y_init = {
        params["y1_init"], params["y2_init"], params["y3_init"], params["y4_init"]
    };


    cop::solver_RungeKutta4<std::valarray<double>, decltype(dglfunc_doppelpendel)*> s_rk4(
        dglfunc_doppelpendel, params["t_init"], y_init
    );

	// write header line
	std::cout<<
        "Time"<<"  "<<
        "\\theta1"<<"  "<<
        "\\theta2"<<"  "<<
        "\\ptilde1"<<"  "<<
        "\\ptilde2"<<"  "<<
        "energy"<<'\n';

    // show parameters
    for (auto p_it = params.begin(); p_it != params.end(); ++p_it)
        std::cout<<"# "<<p_it->first<<" = "<<p_it->second<<'\n';

    // set resolution variable
    std::size_t resolution = params["resolution"], stepcount = 0;

    std::cout<<std::scientific;
    while(s_rk4.getTime() < t_end)
    {
    
        if (! (stepcount % resolution))
        {
            std::cout<<s_rk4.getTime()<<"  "
                <<s_rk4.getY()[0]<<"  "<<s_rk4.getY()[1]<<"  "
                <<s_rk4.getY()[2]<<"  "<<s_rk4.getY()[3]<<"  "
                <<'\n';
        }

        s_rk4.step(timestep);    
        ++stepcount;
    }

    return 0;
}

