#pragma once
#include <utility>
#include <cmath>
#include <vector>
#include <functional>
#include <fstream>

namespace OED_Sover
{    
    template <int arg_num>
    std::array<double, arg_num> Euler(const std::array<std::function<double(std::array<double, arg_num>)>, arg_num-1>& equations,
                                      const std::array<double, arg_num>& init_values, 
                                      double h)
    {
        int equation_num = arg_num-1;
        std::array<double, arg_num> y4;
        for(int i = 0;i<equation_num;++i)
            y4[i+1] = init_values[i+1] + h*equations[i](init_values);
        
        y4[0] = init_values[0] + h;
        return y4;
    }
    
    template <int arg_num>
    std::array<double, arg_num> Classical_Runge_Kutta(const std::array<std::function<double(std::array<double, arg_num>)>, arg_num-1>& equations,
                                                      const std::array<double, arg_num>& init_values, 
                                                      double h)
    {
        int equation_num = arg_num-1;
        double k[equation_num][4];
        std::array<double, arg_num> argument;

        for(int num = 0;num<equation_num;++num)
        {
            k[num][0] = h*equations[num](init_values);
        }
        
        argument[0] = init_values[0] + 0.5*h;
        for(int i = 1;i<arg_num;++i)
            argument[i] = init_values[i] + 0.5*k[i-1][0];
        for(int num = 0;num<equation_num;++num)
        {   
            k[num][1] = h*equations[num](argument);
        }

        for(int i = 1;i<arg_num;++i)
            argument[i] = init_values[i] + 0.5*k[i-1][1];
        for(int num = 0;num<equation_num;++num)
        {   
            k[num][2] = h*equations[num](argument);
        }

        argument[0] = init_values[0] + h;
        for(int i = 1;i<arg_num;++i)
            argument[i] = init_values[i] + k[i-1][2];
        for(int num = 0;num<equation_num;++num)
        {   
            k[num][3] = h*equations[num](argument);
        }

        std::array<double, arg_num> y4;

        for(int i = 0;i<equation_num;++i)
            y4[i+1] = init_values[i+1] + 1/6.0*k[i][0] + 1/3.0*k[i][1] + 1/3.0*k[i][2] + 1/6.0*k[i][3]; 
        
        y4[0] = init_values[0] + h;

        return y4;
    }

    template <int arg_num>
    std::array<double, arg_num> Runge_Kutta_Fehlberg(const std::array<std::function<double(std::array<double, arg_num>)>, arg_num-1>& equations,
                                                     const std::array<double, arg_num>& init_values, 
                                                     double& h, 
                                                     double error_ratio)
    {
        int equation_num = arg_num-1;
        double k[equation_num][6];
        std::array<double, arg_num> argument;

        for(int num = 0;num<equation_num;++num)
        {
            k[num][0] = h*equations[num](init_values);
        }
        
        argument[0] = init_values[0] + 0.25*h;
        for(int i = 1;i<arg_num;++i)
            argument[i] = init_values[i] + 0.25*k[i-1][0];
        for(int num = 0;num<equation_num;++num)
        {   
            k[num][1] = h*equations[num](argument);
        }

        argument[0] = init_values[0] + 3/8.0*h;
        for(int i = 1;i<arg_num;++i)
            argument[i] = init_values[i] + 3/32.0*k[i-1][0] + 9/32.0*k[i-1][1];
        for(int num = 0;num<equation_num;++num)
        {   
            k[num][2] = h*equations[num](argument);
        }

        argument[0] = init_values[0] + 12/13.0*h;
        for(int i = 1;i<arg_num;++i)
            argument[i] = init_values[i] + 1932/2197.0*k[i-1][0] - 7200/2197.0*k[i-1][1] + 7296/2197.0*k[i-1][2];
        for(int num = 0;num<equation_num;++num)
        {   
            k[num][3] = h*equations[num](argument);
        }

        argument[0] = init_values[0] + h;
        for(int i = 1;i<arg_num;++i)
            argument[i] = init_values[i] + 439/216.0*k[i-1][0] - 8.0*k[i-1][1] + 3680/513.0*k[i-1][2] - 845/4104.0*k[i-1][3];
        for(int num = 0;num<equation_num;++num)
        {   
            k[num][4] = h*equations[num](argument);
        }
        

        argument[0] = init_values[0] + 0.5*h;
        for(int i = 1;i<arg_num;++i)
            argument[i] = init_values[i] -8/27.0*k[i-1][0] + 2.0*k[i-1][1] - 3544/2565.0*k[i-1][2] + 1859/4104.0*k[i-1][3] - 11/40.0*k[i-1][4];
        for(int num = 0;num<equation_num;++num)
        {   
            k[num][5] = h*equations[num](argument);
        }

        std::array<double, arg_num> y4;
        std::array<double, arg_num> y5;      

        double temp_min = 1e5;

        for(int i = 0;i<equation_num;++i)
        {
            y4[i+1] = init_values[i+1] + 25/216.0*k[i][0] + 1408/2565.0*k[i][2] + 2197/4104.0*k[i][3] - 0.2*k[i][4]; 

            y5[i+1] = init_values[i+1] + 16/135.0*k[i][0] + 6656/12825.0*k[i][2] + 28561/56430.0*k[i][3] - 9/50.0*k[i][4] + 2/55.0*k[i][5];

            temp_min = std::min(temp_min,(std::abs(init_values[i+1])+std::abs(k[i][0]))/std::abs(y5[i+1]-y4[i+1]));
        }
        
        y4[0] = init_values[0] + h;
        h = 0.9*std::pow(temp_min*error_ratio, 0.2)*h;

        return y4;
    }
}