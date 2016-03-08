# Convecterm
See some nice convection in your terminal:

![Convection in action!](http://i.imgur.com/IjZpNia.png?1 "Convecterm")

Convecterm uses the stokes-flow equations with eulerian advection-diffusion to give you some nice convective behavior in your terminal. It has temperature-dependent viscosity, and (obviously) temperature-dependent densities. It uses finite-differences in a pretty stock-standard way. It does not use any dynamic memory calls, which was a deliberate choice.

It's pretty sensitive to parameter changes, so take it easy when playing with them.

## To run Convecterm:
1. Make sure the CC variable in the Makefile is set to your preferred compiler
2. then just: ```make```
3. ```./convecterm```

## TODO
- Figure out why the dt is so sensitive. Anything much < or > than 1e-4 for the dt makes the model not solve.
