"""
Single Pendulum Simulation.
=================================================

This function simulates a magnetic, base-excited, single pendulum. Both the magnetic and base excitation parts are optional
through setting certain parameters to zero.
"""

def pendulum_simulation(t, IC, params, plotting = True): 
    """This function takes system parameters for a base excited (or not) magnetic (or not) single pendulum
    and simulates a response for desired time frame.
    
    Args:
       t (array):  Time array (1d) strting at 0.
       IC (array): 1d array of length 2 of initial conditions as [initial angle, initial angular velocity].
       params (list): list of float values for parameters. See example below for description of parameters.
       
    Kwargs:
       plotting (bool): Plotting for user interpretation. defaut is False.

    Returns:
       t (array):  Time array (1d) strting at 0.
       th (array):  Array of angles (in degrees) over time t.
       om (array):  Array of anglular velocities (in degrees per second) over time t.

    """
    from scipy.integrate import odeint
    if 'F_phidd' not in globals():
        global F_phidd
        import sympy as sp
        print('Running sympy for equation of motion. May take a few seconds.')
        
        m, l, g, r_cm, I_cm, A, w, mu_v, mu_q, mu_c, q, d, mu = sp.symbols("m \ell g r_{cm} I_{cm} A \omega \mu_v \mu_q \mu_c q d \mu") #declare constants
        t_s = sp.symbols("t") #declare time variable
        th = sp.Function(r'\theta')(t_s) #declare time dependent variables
        
        r = sp.sqrt((l)**2 +(d+l)**2 - 2*l*(l+d)*sp.cos(th))
        a = 2*np.pi - np.pi/2
        b = (np.pi/2) - th
        phi = np.pi/2 - sp.asin((l/r)*sp.sin(th))
        Fr = (3*mu*q**2/(4*np.pi*r**4))*(2*sp.cos(phi-a)*sp.cos(phi-b) - sp.sin(phi-a)*sp.sin(phi-b))
        Fphi =  (3*mu*q**2/(4*np.pi*r**4))*(sp.sin(2*phi-a-b))
        
        tau_m = l*Fr*sp.cos(phi-th) - l*Fphi*sp.sin(phi-th)
        tau_v = mu_v*th.diff(t_s)
        tau_q = mu_q*(th.diff(t_s)**2)*th.diff(t_s)/(np.abs(th.diff(t_s))+10E-6)
        tau_c = mu_c*th.diff(t_s)/(np.abs(th.diff(t_s))+10E-6)
        
        V = -m*g*r_cm*sp.cos(th)
        vx = r_cm*th.diff(t_s)*sp.cos(th) + A*w*sp.cos(w*t_s)
        vy = r_cm*th.diff(t_s)*sp.sin(th)
        T = (1/2)*I_cm*th.diff(t_s)**2 + (1/2)*m*(vx**2 + vy**2)
        R = tau_v + tau_q + tau_c + tau_m
        
        L = T - V
        R_symb = sp.symbols("R_s")
        EOM = (L.diff(th.diff(t_s))).diff(t_s) - L.diff(th) + R_symb #lagranges equation applied
    
        #first solve both EOM for th1dd and th2dd
        thdd = sp.solve(EOM, th.diff(t_s).diff(t_s))[0]
    
        #we first need to change to th_1 and om_1 symbols and not functions to apply lambdify.
        ph, phi_dot = sp.symbols(r"\phi \dot{\phi}")
        phidd = thdd.subs([(R_symb,R)])
        phidd = phidd.subs([(th.diff(t_s),phi_dot)])
        phidd = phidd.subs([(th,ph)])
        
        #lambdified functions
        F_phidd = sp.lambdify([(t_s, ph, phi_dot), (m, l, g, r_cm, I_cm, A, w, mu_v, mu_q, mu_c, q, d, mu)], phidd)
        print('sympy is complete!')
        
        
    
    p = params
    w0 = IC
    def vectorfield(w, t, p):
        ph, phi_dot = w
        M, l, g, r_cm, I_cm, A, w, mu_v, mu_q, mu_c, m, d, mu = p
        f = [phi_dot,
             F_phidd((t, ph, phi_dot), (M, l, g, r_cm, I_cm, A, w, mu_v, mu_q, mu_c, m, d, mu))]
        return f
    
    sol = odeint(vectorfield, w0, t, args=(p,))
    rad2deg = 180/np.pi
    th, om = sol.T[0]*rad2deg, sol.T[1]*rad2deg
    t = t
    
    if plotting == True:
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec

        gs = gridspec.GridSpec(2,4) 
        plt.figure(0, figsize=(10,4.5)) 
        
        plt.subplot(gs[0, 0:2])
        plt.ylabel(r'$\theta$ ($^\circ$)')
        plt.plot(t, th, 'b')
        
        plt.subplot(gs[1, 0:2])
        plt.ylabel(r'$\dot{\theta}$ ($^\circ$/s)')
        plt.xlabel('time (s)')
        plt.plot(t, om, 'g')
        
        plt.subplot(gs[0:2, 2:4])
        plt.xlabel(r'$\theta$ ($^\circ$)')
        plt.ylabel(r'$\dot{\theta}$ ($^\circ$/s)')
        plt.plot(th, om, 'k')
        
        plt.subplots_adjust(wspace=0.7)
        plt.savefig('C:\\Users\\myersau3.EGR\\Desktop\\python_png\\mag_unknown_sing_pend_simul.png', bbox_inches='tight',dpi = 400)
        plt.show()   
    
    return t, th, om
   
# In[ ]: 

if __name__ == "__main__": #only runs if executed through this script
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
    
        
        