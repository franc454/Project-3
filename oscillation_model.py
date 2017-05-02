# -*- coding: utf-8 -*-
"""
Project 3

This is a model for harmonic oscillation; both simple and damped.
Euler's method is used here.  Euler's method begins to fail when dt > 1 and the
location of the oscillating mass will continue to grow to positive and negative
infinity.  A better method would be a Runge-Kutta method.

Author: Keenen
"""

from matplotlib import pyplot as plt
import numpy as np
import math as mt
from matplotlib import animation

k = 1  #some constant
m = 1  #mass
c = 0.1  #viscous damping coefficient

#Damping ratio
#Underdamped: d_ratio < 1
#Critically damed: d_ratio = 1
#Overdamped: d_ratio > 1
d_ratio = c/(2*mt.sqrt(m*k))

#undamped angular frequency of the oscillator
ang_freq = mt.sqrt(k/m)

#Be careful with the time interval (dt).  Euler method begins to fail when dt > 1.
dt = 0.5  #time interval
t = 0  #current time
total_t = 1000  #total time

#make empty lists that are the length of total time divided by the dt.  This 
#is so that we can place the next location of the mass at each time step.  We
#will later plot these locations to visualize the oscillation.
x_prev = np.empty([int(total_t/dt)+1])
x_now= np.empty([int(total_t/dt)+1])

x_prev_damp = np.empty([int(total_t/dt)+1])
x_now_damp = np.empty([int(total_t/dt)+1])

#create a list the same length as above lists but filled with some constant so 
#the mass oscillates along one x-value on a plot
y = np.empty([int(total_t/dt)+1])
y.fill(10)

#setting initial timesteps so that way motion has been initiated
x_prev[0] = 100
x_now[0] = 99

x_prev_damp[0] = 100
x_now_damp[0] = 99

#value to step through lists of mass locations sowe can assign locations.
i = 0

#x_next, x_prev and x_now are used in the Euler equation for the simple 
#harmonic oscillator. x_next_damp, x_prev_damp and x_now_damp are used for the
#damped harmonic oscillator.  These equations calculate the location of the mass
#at the next time interval based on the current and previous locations.  
while t < total_t:
    x_next = ((-k/m)*x_now[i]*(dt)**2)+(2*x_now[i])-x_prev[i]
    x_prev[i+1] = x_now[i]
    x_now[i+1] = x_next
    
    x_next_damp = ((2*d_ratio*dt*ang_freq*x_now_damp[i])-((dt**2)*(ang_freq**2)*\
            x_now_damp[i])-x_prev_damp[i]+(2*x_now_damp[i]))/(2*d_ratio*dt*ang_freq + 1)
    x_prev_damp[i+1] = x_now_damp[i]
    x_now_damp[i+1] = x_next_damp
    
    t += dt
    i += 1
    
#create figure with 2 plots so we can see both the simple and damped oscillations.  
fig = plt.figure(figsize = (10,10))

#simple harmonic oscillator plot
graph = fig.add_subplot(211)
graph.set_title('Simple Harmonic Oscillator')
point, = graph.plot([x_now[0]], [y[0]], 'ro')
graph.grid()

#damped harmonic oscillator plot
graph_damp = fig.add_subplot(212)
graph_damp.set_title('Damped Harmonic Oscillator')
point_damp, = graph_damp.plot([x_now_damp[0]], [y[0]], 'ro')
graph_damp.axes.get_yaxis().set_visible(False)
graph_damp.grid()

#set x-axis limits so the mass is within a reasonable view     
graph.set_xlim([-x_prev[0]*1.5, x_prev[0]*1.5])  
graph_damp.set_xlim([-x_prev[0]*1.5, x_prev[0]*1.5])  
graph.axes.get_yaxis().set_visible(False)
#returns a new point to plot on the same plot.  This removes the previously 
#plotted points so that the plot is animated and we can actually watch the 
#mass oscillate.  We're essentially taking the plot and swtiching the plotted 
#point value to the next value in the x_now or x_now_damped list.
def oscillate(n, x_now, x_prev, point):
    point.set_data(np.array([x_now[n], y[n]]))
    return point

def damped_oscillate(n, x_now_damp, x_prev_damp, point_damp):
    point_damp.set_data(np.array([x_now_damp[n], y[n]]))
    return point_damp

#calling each definition to switch plotted point for so many frames and create
#the animation we want.
ani = animation.FuncAnimation(fig, oscillate, frames = int(total_t/dt), \
                        fargs=(x_now, x_prev, point))
ani_damp = animation.FuncAnimation(fig, damped_oscillate, frames = int(total_t/dt), \
                        fargs=(x_now_damp, x_prev_damp, point_damp))

plt.show()