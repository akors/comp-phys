#include <string>
#include <map>
#include <valarray>

#include <iostream>
#include <fstream>
#include <sstream>


#define _USE_MATH_DEFINES
#include <cmath>

// initialize with default paramters
static std::map<std::string, double> params {
        {"stepsize", 0.001},
        {"t_init", 0.0}, // initial conditions
        {"y_init", 1.0},
        {"t_end", 3},
        {"theta0", 0},
        {"dotheta0", 2*M_PI},
        {"k_exp", 1.0}
};


void readParams(std::map<std::string, double>& parms, std::ifstream& parmfile)
{
    std::size_t linecount = 1;
    std::string line;
    std::string n; double v;

    while(!parmfile.eof())
    {
        // read line
        getline(parmfile, line);

        // skip empty lines
        if(line.empty())
            continue;

        // read parameters
        std::stringstream ss(line);
        ss>>n>>v;

        // bail out on errors, store on success
        if (!ss)
        {
            std::clog<<"Error reading parameter in line "<<linecount<<'\n';
            continue;
        }
        else
            parms[n] = v;

        ++linecount;
    }
}

template <typename T, typename F> 
struct solver_euler
{
    double t_cur;
    T y_cur; // current t/y values
    F f; // functor depending on t,y

    // initializing constructor
    solver_euler(const F& func, double t_init, T y_init)
        : f(func), t_cur(t_init), y_cur(y_init)
    {     
    }

    // no default-construction!
    solver_euler() = delete;

    T step(double timestep)
    {
        y_cur = y_cur + f(t_cur, y_cur) * timestep;
        t_cur += timestep;
        return y_cur;
    }
};

template <typename T, typename F> 
struct solver_RungeKutta4
{
    double t_cur;
    T y_cur; // current t/y values
    F f; // functor depending on t,y

    // initializing constructor
    solver_RungeKutta4(const F& func, double t_init, T y_init)
        : f(func), t_cur(t_init), y_cur(y_init)
    {     
    }

    // no default-construction!
    solver_RungeKutta4() = delete;

    T step(double timestep)
    {
        T k1 = f(t_cur             , y_cur                  );
        T k2 = f(t_cur + timestep/2, y_cur + timestep/2 * k1);
        T k3 = f(t_cur + timestep/2, y_cur + timestep/2 * k2);
        T k4 = f(t_cur + timestep  , y_cur + timestep   * k3);

        y_cur = y_cur + timestep/6*(k1 + 2.0*k2 + 2.0*k3 + k4);

        t_cur += timestep;
        return y_cur;
    }
};

#if 0
double dglfunc_exponential(double /* t */, double y)
{
    return params["k_exp"]*y;
}
#endif

std::valarray<double> dglfunc_exponential(double /* t */, std::valarray<double> y)
{
    return params["k_exp"] * y;
}

int main()
{
    std::ifstream parmfile("02_pendulum.prm");
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


    // show parameters
    for (auto p_it = params.begin(); p_it != params.end(); ++p_it)
        std::cout<<"# "<<p_it->first<<" = "<<p_it->second<<'\n';

    auto y_init = std::valarray<double>{params["y_init"]};

    solver_euler<std::valarray<double>, decltype(dglfunc_exponential)*> s_euler(
        dglfunc_exponential, params["t_init"], y_init
    );
    solver_RungeKutta4<std::valarray<double>, decltype(dglfunc_exponential)*> s_rk4(
        dglfunc_exponential, params["t_init"], y_init
    );

    std::cout<<"# Writing data for in ["<<params["t_init"]<<", "<< params["t_end"]<<") with delta = "<<params["stepsize"]<<'\n';
    std::cout<<std::scientific;

    double t = params["t_init"], timestep = params["stepsize"];
    while(t < params["t_end"])
    {
        std::cout<<
                t<<"  "<<
                params["y_init"] * std::exp(params["k_exp"]*t)<<"  "<<
                s_euler.step(timestep)[0]<<"  "<<
                s_rk4.step(timestep)[0]<<
                '\n';
        t += timestep;
    }

    return 0;
}

