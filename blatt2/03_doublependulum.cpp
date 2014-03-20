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
        {"y2_init", 2},
        {"omega", 1}
};


std::valarray<double> dglfunc_doppelpendel(double /* t */, const std::valarray<double>& y)
{
    double d_theta = y[0] - y[1];

    return {
        y[1],
        - params["omega"]*params["omega"] * std::sin(y[0])
    };
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
        readParams(params, parmfile);
        std::clog<<"Parameters loaded.\n";
    }




	// initialize solver
    double timestep = params["stepsize"], t_end = params["t_end"];
    std::valarray<double> y_init = {params["y1_init"], params["y2_init"]};


    solver_RungeKutta4<std::valarray<double>, decltype(dglfunc_doppelpendel)*> s_rk4(
        dglfunc_doppelpendel, params["t_init"], y_init
    );

	// write header line
	std::cout<<
        "Time"<<"  "<<
        "\\theta"<<"  "<<
        "\\dottheta"<<"  "<<
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
                <<s_rk4.getY()[0]<<"  "<<s_rk4.getY()[1]<<"  "<<energy(s_rk4.getY())
                <<'\n';
        }

        s_rk4.step(timestep);    
        ++stepcount;
    }

    return 0;
}

