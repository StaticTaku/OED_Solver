# OED_Solver for C++

## How to use

This program is header only.
You can just use all of integral solver scheme by just including.

```c++
    #include <OED_Solver.hpp>
```

All of N-order OED return to N simultaneous OED equations,so this program only focus on soloving N simultaneous OED equations.

Let the independent variable be $t$,and let the N simultaneous OED equations like this below.

$$
\begin{eqnarray*}
\frac{dy_1(t)}{dt} &=& f_1(t,y_1,y_2,\cdot\cdot\cdot,y_N) \\
\frac{dy_2(t)}{dt} &=& f_2(t,y_1,y_2,\cdot\cdot\cdot,y_N) \\
&\cdot& \\
&\cdot& \\
&\cdot& \\
\frac{dy_N(t)}{dt} &=& f_N(t,y_1,y_2,\cdot\cdot\cdot,y_N)
\end{eqnarray*}
$$

In this program, you can define N simultaneous OED equations like this below

```c++
    const int arg_num = N+1 //N is number of equations, and the number of independent variable is 1 so arg num is
    auto f1 = [](std::array<double, arg_num> args) { return ;};
    auto f2 = [](std::array<double, arg_num> args) { return ;};
    .
    .
    .
    auto fN = [](std::array<double, arg_num> args) { return ;};
```

Next, you have to set initial value of $y_1,y_2,\cdot\cdot\cdot,y_N$ at $t$

In this program, you can define initial value at $t$ like this below

```c++
    std::array<double, arg_num> init_values = {t,y1,y2, ... ,yN};
```

Let $dt>0$ be step size for solving OED

Finally, you can solve OED like this

```c++
    double dt = 1e^-3;
    auto ans = OED_Sover::Euler<arg_num>({f1,f2,...,fN},init_values,dt);
```

Content of ans is like this

```c++
    ans = {t+dt,y1(t+dt),y2(t+dt),...,yN(t+dt)}
```

If you want to solve OED equations until $t = goal$,you can do like this

```c++
    std::array<double, arg_num> init_values = {t,y1,y2, ... ,yN}
    while(t < goal)
    {
        init_values = OED_Sover::Euler<arg_num>({f1,f2,...,fN},init_values,dt);
        t           = init_values[0];
    }
```

you can define array of {f1,f2,...,fN} like this below

```c++
    std::array<std::function<double(std::array<double, arg_num>)>, arg_num-1> equations = {f1,f2,...,fN}
```

## Example

Lets solve this

$$
\begin{eqnarray*}
\frac{dx}{dt} &=& x+2y \\
\frac{dy}{dt} &=& 2x+y \\
\end{eqnarray*}
$$

In c++,

```c++
    const int arg_num = 3;

    auto func1 = [](std::array<double, arg_num> args) { return args[1] + 2*args[2];};
    auto func2 = [](std::array<double, arg_num> args) { return 2*args[1] + args[2];};
```

Let initial values be like this below

```c++
    std::array<double, arg_num> init_values = {0,1,2};
```

and, lets solve this problem in c++ !!

```c++
    while(t < 2)
    {
        init_values = OED_Solver::Euler<arg_num>({func1,func2},init_values,dt);
        t           = init_values[0];
    }
```

OED_Solver::Runge_Kutta_Fehlberg is the scheme that you can decide accuracy of answer
you can use it like this below

```c++
    real error_ratio = 1e-5;
    real ans = OED_Solver::Runge_Kutta_Fehlberg<arg_num>({f1,f2,...,fN},init_values,dt,error_ratio);
```

You can define error_ratio of ans comparing exact ans. dt is going to be updated into another value.
