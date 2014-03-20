#ifndef COPLIB_HPP_INCLUDED
#define COPLIP_HPP_INCLUDED

#include <map>
#include <string>
#include <sstream>


namespace cop {

void readParams(std::map<std::string, double>& parms, std::ifstream& parmfile);

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
        T k1 = f(t_cur               , y_cur                    );
        T k2 = f(t_cur + timestep/2.0, y_cur + timestep/2.0 * k1);
        T k3 = f(t_cur + timestep/2.0, y_cur + timestep/2.0 * k2);
        T k4 = f(t_cur + timestep    , y_cur + timestep     * k3);

        y_cur = y_cur + timestep/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);

        t_cur += timestep;
        return y_cur;
    }
};

} // namespace cop

#endif // ifndef COPLIB_HPP_INCLUDED

