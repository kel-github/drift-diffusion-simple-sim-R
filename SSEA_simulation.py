import numpy as np
from bokeh.io import output_notebook, show
from bokeh.plotting import figure

def diffusion(dv: float, 
              v: float,
              a: float,
              s: float, 
              t: float,
              dt: float, 
              t_er: float) -> tuple:
    '''
    Simulate one diffusion time trace until a decision criterion is reached.
    This simulation uses Euler’s method to generate sample paths of the diffusion process.
    This is the loop that controls the evidence accumulation process. 
    At each time step, the current evidence total (dv) is compared against the value of the two absorbing boundaries (a and 0). If the evidence total is sitting between these to values the process has not yet terminated (i.e., neither boundary has been reached) and so more information is sampled. 
    Each sample of information changes the evidence total by an amount that is fixed across all time steps (v * dt; or the mean Drift Rate multiplied by the size of the time step, here set to 1 ms) plus a (scaled) random draw from a normal distribution with mean 0 and standard deviation equal to s. 
    If the process has reached/exceeded one of the absorbing boundaries, the number of time steps required to reach that threshold is recorded along with the absorbing boundary that was reached.
    Params:
        - dv: decision variable
        - v: drift rate
        - a: boundary separation
        - s: within-trial variability
        - t: time in seconds
        - dt: time step in seconds
    Returns:
        - tuple containing the response and the time it took to make a decision
    '''
    while dv > 0 and dv < a:
        dv += v * dt + s * np.random.randn() * np.sqrt(dt)
        t += dt
    rsp = 0
    if dv > 0: 
        rsp = 1
    RT = t + t_er
    return (rsp, RT)

def simulateDiffusion(n_trials = 1000,
                      v = .2,
                      a = .04,
                      t_er = .2,
                      bias = .5,
                      s = .1,
                      dt = .001):
    '''
    Quick and dirty method of generating a bunch of simulated choice-RT
    Params:
        - n_trials: how many trials should be simulated
        - v: drift rate
        - a: boundary separation
        - t_er: .2
        - bias: starting point / response bias
        - s: diffusion coefficient / within-trial standard deviation
        - dt: time step in seconds

    The s parameters is traditionally fixed to 0.1 (following Ratcliff, 1978), but others set s = 1. Either is fine, as this parameter “scales” the other parameters in the model (i.e., diffusion model parameters are only defined on a ratio scale, and so you can get the same predictions from different combinations of parameter values - if one parameter were to be mulitplied by a constant, the same predictions would obtain if all other parameters were multiplied (or divided) by the same value).

    ''' 
    z = a * bias
    return np.array([diffusion(z, v, a, s, 0, dt, t_er) for trial in np.arange(n_trials)])
 
response, RT = simulateDiffusion().T

p = figure(plot_height=300)
histA, edgesA = np.histogram(RT[response == 1], density=True, bins=20)
p.quad(top=histA, bottom=0, left=edgesA[:-1], right=edgesA[1:], 
       alpha=0.4, fill_color = 'blue', line_color = 'gray', legend = 'Correct')
histB, edgesB = np.histogram(RT[response == 0], density = True, bins = 20)
p.quad(top = histB, bottom = 0, left = edgesB[:-1], right = edgesB[1:], 
       alpha = 0.4, fill_color = 'orange', line_color = 'gray', legend = 'Error')

output_notebook()
show(p)
