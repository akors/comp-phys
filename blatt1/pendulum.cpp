#include <string>
#include <map>

#include <iostream>
#include <fstream>
#include <sstream>



#define _USE_MATH_DEFINES
#include <cmath>

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

int main()
{
    // initialize with default paramters
    std::map<std::string, double> params {
            {"stepsize", 0.001},
            {"theta0", 0},
            {"dotheta0", 2*M_PI}
    };

    std::ifstream parmfile("params.txt");
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


    return 0;
}

