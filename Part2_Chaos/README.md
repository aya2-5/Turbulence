# EPFL ME467: Turbulence
### Time-marching the KSE
The following lines show how to utilize the function `KSE_integrate` to advance the KSE in time.
```
close all; clear; clc

addpath('../functions/')
L = 38.6;
N = 64;
symm = true;
T = 750;
dt = 0.1;

[x,~] = domain(L,N);
u0 = sin(2*pi*x/L);

v0 = field2vector(u0,N,symm);
[vT,~] = KSE_integrate(v0,T,dt,0,L,N,symm);

uT = vector2field(vT,N,symm);

figure
  plot(x,uT,'LineWidth',2)
  hold on; grid on
  plot(x,u0,'LineWidth',2)
  xlabel('x'); ylabel('u')
  legend('u(T)','u(0)')
```
This program can be broken down into the following steps:
- We first tell MATLAB to look for functions in the relative folder path `../functions/`.
- Then, we define the domain length `L = 38.6`, number of grid points `N = 64`, imposed center symmetry `symm = true`, integration time `T = 750` and time step size `dt = 0.1`.
- To define the initial condition, we first construct the grid points using `[x,~] = domain(L,N)`. Then, the initial condition is defined using the sine function `u0 = sin(2*pi*x/L)`. Note that this initial condition is consistent with the imposed center symmetry.
- We construct the state vector associated with `u0` by calling `v0 = field2vector(u0,N,symm)`. Finally, we call `[vT,~] = KSE_integrate(v0,T,dt,0,L,N,symm)` where `vT` is the state vector of the KSE obtained by advancing `v0` for time `T`.
- Since `vT` is a state vector and not the vector of grid values, we call `uT = vector2field(vT,N,symm)` to transform `vT` to its corresponding physical field `uT` and, then, plot `u0` and `uT` in one figure.

 **Note 1**: By setting the 4th argument of `KSE_integrate` to `0`, no intermediate snapshot is stored and `vT` is a column vector. In order to store intermediate snapshots, you can call `[vt,t] = KSE_integrate(v0,T,dt,dt_store,L,N,symm)`. As a result, every `dt_store` the snapshot is appended to `vt` as a new column. The vector `t` contains the time instance associated with each column of `vt`.

### 3D projection of the state space
The following lines show how to use the function `projection` to visualize a 3D projection of the chaotic trajectory from the previous example.
```
close all; clear; clc

addpath('../functions/')
L = 38.6;
N = 64;
symm = true;
T = 750;
dt = 0.1;
dt_store = 1;

[x,~] = domain(L,N);
u0 = sin(2*pi*x/L);

v0 = field2vector(u0,N,symm);
[vv,tt] = KSE_integrate(v0,T,dt,dt_store,L,N,symm);

[E,P,D] = projection(vv,N,L,symm);

figure
  plot3(E,P,D,'Color',[0,0,0,0.25])
  xlabel('E'); ylabel('P'); zlabel('D'); grid on
```
Here `[0,0,0,0.25]` specifies red, green, blue and the opacity of the line. We have chosen a semi-transparent line so that regions of the state space that are visited more often appear darker than the regions that are rarely visited.

### Search for a periodic orbit or an equilibrium
The following program is an example of how to utilize the function `search4PO` to compute a periodic orbit from a guess. Here, the same sine function as in the previous examples is used, and a period of `T_guess=50` is chosen for demonstration.
```
close all; clear; clc

addpath('../functions/')
L = 38.6;
N = 64;
symm = true;
dt = 0.1;

[x,~] = domain(L,N);
u0 = sin(2*pi*x/L);

v_guess = field2vector(u0,N,symm);
T_guess = 50;

[v_best,T_best,flag] = search4PO(v_guess,T_guess,dt,L,N,symm); 
```
This crude guess does not converge, and the following will be displayed:
```
                                         Norm of      First-order   Trust-region
 Iteration  Func-count     f(x)          step         optimality    radius
     0         23         4146.85                      2.05e+03               1
     1         46         3586.44              1            600               1
   ...        ...             ...            ...            ...             ...
    75       1594         2927.01      0.0149012           63.7          0.0149

Solver stopped prematurely.

fsolve stopped because it exceeded the iteration limit,
options.MaxIterations = 7.500000e+01.
```
In a successful search, the search flag takes the value of 1. Here, however, `flag = 0`. If the search is successful, `v_best` is the converged state vector, and `T_best` is the converged period of the periodic orbit. In that case, if you advance `v_best` as the initial condition for a time interval of `T_best` (see the first example above), the trajectory closes on itself. Note that `v_guess` and `v_best` are both state vectors.

In order to search for an equilibrium solution, we use `[v_best,flag] = search4EQ(v_guess,T_eqb,dt,L,N,symm)`, where `T_eqb` is an arbitrary integration time.
