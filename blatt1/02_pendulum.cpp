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
        {"k_exp", 0.33}
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
    T t_cur, y_cur; // current t/y values
    F f; // functor depending on t,y

    // initializing constructor
    solver_euler(const F& func, T t_init, T y_init)
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

double dglfunc_exponential(double /* t */, double y)
{
    return params["k_exp"]*y;
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

    solver_euler<double, decltype(dglfunc_exponential)*> s_euler(
        dglfunc_exponential, params["t_init"], params["y_init"]
    );

    std::cout<<"# Writing data for in ["<<params["t_init"]<<", "<< params["t_end"]<<") with delta = "<<params["stepsize"]<<'\n';
    std::cout<<std::scientific;

    double t = params["t_init"], timestep = params["stepsize"];
    while(t < params["t_end"])
    {
        std::cout<<t<<"  "<<s_euler.step(timestep)<<'\n';
        t += timestep;
    }

    return 0;
}

