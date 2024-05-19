# linadv_ser.jl
* **Linear advection equation** is a **hyperbolic partial differential equation** describes the motion of a scalar quantity `u` as it is advected by a known velocity field. 

* **`linadv_ser.jl` uses Lax-Wendroff method for updating the solution at each time step until the final time.**   

* This folder contain code for solving **linear advection equation** in 1D in a *serial* programming approach.  

## Usage
* Use the following to compile and `linadv_ser.jl` file. This operation will create 2 files in the folder whose path you can set in code of `linadv_ser.jl`  
*   ```
    julia path/to/file/linadv_ser.jl
    ```  
* Two files: `exact_sol.txt` and `num_sol.txt` will get created, contains *exact solution* and *numerial solution* respectively.

* Run plotting script to visualize these solution:
*   ```
    julia plot.jl
    ```

## Equation
* Form of the equation is:  `u_t + a*u_x = 0`  
where-  
`u_t` is partial derivative with respect to `t`  
`u_x` is partial derivative with respect to `x`  
`a` is advection velocity  

* Computational domain: `[0,1]`

## Initial Conditions
* Following initial conditions have been used here: `u(x) = sin(2πx)`

## Functions
This section explains the `functions` involved in the `linadv_ser.jl` code.  
* `initial_u!(param::NamedTuple, x::mesh, u::Vector{Float64})`
    * This function generates initial solution vector based on the initial conditions provided.

* `exact_solution!(param::NamedTuple, x::mesh, exact_sol::Vector{Float64})`  
    * This function generates exact solution vector of the problem (if expression is known to the user) to find the error.

*  `update_lw!(u::Vector{Float64}, unew::Vector{Float64}, sigma::Float64)`
    * This function is implements `Lax-wendroff` method and updates the solution at different time levels.

    * Mathematical description:

    $$u_i^{j+1} = u_i^{j} - \frac{\sigma}{2} (u_{i+1}^{j} - u_{i-1}^{j}) + \frac{\sigma}{2}^2 (u_{i-1}^j - 2 u_i^j + u_{i+1}^j)$$  

    where: $\Delta t$ = time stepping; $\Delta x$ = grid spacing (divide by `N-1` not `N` to avoid problem with boundary conditions)

    $$\sigma == sigma == a \frac{\Delta t}{\Delta x}$$  

* `error_cal(param::NamedTuple, exact_sol::Vector{Float64}, u::Vector{Float64})`
    * This function calculates and return absolute error in numerical solution.

* `solver(param::NamedTuple)`  
    * This function contains the main functionality of the code and perform the calling to the functions that have been described above.

    * This function will perform a input operation of `cfl` from the user.

    * It prints out the `Iteration` and `Error` to the `stdout`.  

    * It outputs the numerical and exact solution to `.txt` files that will be used for plotting the solution in `plot.jl`  

    * `Important Point`: Inside `while` loop of this function, there is a `crucial block`. Why is it so?  
        * While looping from `t=0` to `t=t_final`, the looping index `j` can go beyond `t_final` ($j + dt > t$) due to the so value of `dt`.
        * Implemented `if` condition is needed to check the above and if `if` is `true` then `dt` needs to be changed to $dt = t - j$ and hence `sigma` needs to be changed accordingly.
## Plotting
* See `Usage` section for plotting of numerical solution along with the exact solution.  
* Graph will be in the `.png` file.