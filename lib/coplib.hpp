#ifndef COPLIB_HPP_INCLUDED
#define COPLIP_HPP_INCLUDED

#include <map>
#include <string>
#include <sstream>


namespace cop {

typedef double real_t;

void readParams(std::map<std::string, real_t>& parms, std::ifstream& parmfile);

template <typename T, typename F> 
struct solver_euler
{
    real_t t_cur;
    T y_cur; // current t/y values
    F f; // functor depending on t,y

    // initializing constructor
    solver_euler(const F& func, real_t t_init, T y_init)
        : f(func), t_cur(t_init), y_cur(y_init)
    {     
    }

    // no default-construction!
    solver_euler() = delete;

    T step(real_t timestep)
    {
        y_cur = y_cur + f(t_cur, y_cur) * timestep;
        t_cur += timestep;
        return y_cur;
    }
};


template <typename T, typename F> 
struct solver_RungeKutta4
{


    // initializing constructor
    solver_RungeKutta4(const F& func, real_t t_init, T y_init)
        : f(func), t_cur(t_init), y_cur(y_init)
    {     
    }

    // no default-construction!
    solver_RungeKutta4() = delete;

    void step(real_t timestep)
    {
        k1 = f(t_cur               , y_cur                    );
        k2 = f(t_cur + timestep/2.0, y_cur + timestep/2.0 * k1);
        k3 = f(t_cur + timestep/2.0, y_cur + timestep/2.0 * k2);
        k4 = f(t_cur + timestep    , y_cur + timestep     * k3);

        y_cur += timestep/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);

        t_cur += timestep;
    }
    
    const T& getY() const
    { return y_cur; }
    
    const real_t& getTime() const
    { return t_cur; }

private:
    // store these locally, so we don't have to create the k's on every step
    T k1, k2, k3, k4;

    real_t t_cur;
    T y_cur; // current t/y values
    F f; // functor depending on t,y
};

} // namespace cop

#endif // ifndef COPLIB_HPP_INCLUDED

