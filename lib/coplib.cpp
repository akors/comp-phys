#include "coplib.hpp"

#include <iostream>
#include <fstream>

using namespace cop;

void cop::readParams(std::map<std::string, double>& parms, std::ifstream& parmfile)
{
    std::size_t linecount = 0;
    std::string line;
    std::string n; double v;

    while(!parmfile.eof())
    {
        ++linecount;

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
    }
}

