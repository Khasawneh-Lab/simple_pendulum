# simple_pendulum

This repo contains (or links to) the simulation Python code, CAD, model, BOM, mechanical drawings, equations, and parameters for a manufactured single pendulum. A CAD rendering of the model is provided below.

<p align="center">
  <img src="https://github.com/Khasawneh-Lab/simple_pendulum/blob/master/figures/single_pendulum_fig.png">
</p>

## Download

The CAD will soon be available through GrabCad.

## Usage and Documentation

Full documentation of the pendulum is provided [here](https://github.com/Khasawneh-Lab/simple_pendulum/blob/master/simple_pendulum_documentation.pdf).

This material is based upon work supported by the National Science Foundation under grant numbers CMMI-1759823 and DMS1759824. 

Disclaimer: Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.

## Simulation Documentation

    #import packages
    import numpy as np
    
    #----------------------------User Defined Parameters------------------------------------------
    # forcing parameters
    w, A = 3*np.pi, 0.0097 
    # description: base frequency, base amplitude
    # units: rad/s, m
    
    #magnetic parameters
    m, mu, d = 0.85, 1.257E-6, 0.036 
    # description: magnetic dipole magnitude, free space permeability, distance
    # units: C*m, ?, m
    
    #damping parameters
    mu_v, mu_q, mu_c = 0.00012, 0.000009, 0.0022 
    # description: viscous, quadratic, and coulomb damping constants
    # units: Ns/m, Ns^2/m^2, N
    
    #pendulum parameters
    M, l, r_cm, I_cm, g = 0.1038, 0.208, 0.142, 0.00071, 9.81 
    # description: mass, length, distance to center of mass, Inertia about center of mass, gravity
    # units: kg, m, m, kg*m^2, m/s^2
    
    #Pack parameters into list
    params = [M, l, g, r_cm, I_cm, A, w, mu_v, mu_q, mu_c, m, d, mu]
    
    #Initial conditions
    th_0, om_0 = 2.0, 0.0 # rad, rad/s
    IC = [th_0, om_0]
    
    #ODE solve time
    stoptime, fs = 40, 100
    t = np.linspace(0,stoptime,stoptime*fs)
    #------------------------------------------------------------------------------------------
    
    
    #run simulation functions to get time array, angle (th), and angular veloctiy (om)
    t, th, om = pendulum_simulation(t, IC, params, plotting = True)


<p align="center">
  <img src="https://github.com/Khasawneh-Lab/simple_pendulum/blob/master/figures/simulation_fig.png">
</p>

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
