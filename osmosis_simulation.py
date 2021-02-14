# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 14:15:36 2021

@author: Vivian Imbriotis
"""

import numpy as np
from numpy.linalg import norm as mod
from numpy.linalg import inv
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import seaborn as sns
from osmosisbackend import collision_1d
from itertools import chain

sns.set_style("dark")

RNG = np.random.default_rng()
LOWER = np.array((0,0))
UPPER = np.array((1,1))
UPPER_SOLUTE = np.array((0.5,1))

def create_particle_at_temp(temp,lower_bounds=LOWER,upper_bounds=UPPER):
    pos = RNG.random(2) * (upper_bounds - lower_bounds) - lower_bounds
    vel = RNG.random(2) * 2 - 1
    #modalize
    vel /= (vel**2).sum()**0.5
    #temp
    vel *= temp
    return Particle(pos,vel,xmin=lower_bounds[0],xmax=upper_bounds[0],
                    ymin=lower_bounds[1],ymax= upper_bounds[1])

def create_list_of_particles_at_temp(N,temp,lower_bounds=LOWER,upper_bounds=UPPER):
    return [create_particle_at_temp(temp,lower_bounds,upper_bounds) for _ in range(N)]


##def collision_1d(u1, u2, m1, m2):
##    
##    v1 = (m1-m2)/(m1+m2)*u1 + (2*m2)/(m1+m2)*u2
##    v2 = 2*m1/(m1+m2)*u1 + (m2-m1)/(m1+m2)*u2
##    return (v1,v2)
    

#Eventually we're going to need different kinds of particles that have a mask
#for which membranes with which to collide.
class Particle:
    resistivity = 0.005
    repulsion = 0.05
    m = 1 #mass
    def __init__(self,pos,vel,xmin=0,xmax=1,ymin=0,ymax=1):
        self.pos = pos
        self.vel = vel
        self.collision_handled = 0
        self.xmin, self.ymin, self.xmax,self.ymax = xmin,ymin,xmax,ymax
    def update(self,ls_of_particles):
        if not self.collision_handled:
            self.check_collision(ls_of_particles)
        else: 
            self.collision_handled -=1
        if self.xmin>self.pos[0] or self.xmax<self.pos[0]: self.vel[0] *= -1
        if self.ymin>self.pos[1] or self.ymax<self.pos[1]: self.vel[1] *= -1

        self.pos += self.resistivity*self.vel
    def check_collision(self,ls_of_particles):
        for p in ls_of_particles:
            if p is not self and mod(p.pos - self.pos) < self.repulsion:
                delta_pos = p.pos - self.pos
                #We need to seperate the tangent and normal components
                # of the collision. The tangent component of both
                # particles is unchanged, and the normal component is
                # updated as per a 1-dimentional collision (i.e. in accord
                # with cons of energy and momentum!)
                norm_unit_v = delta_pos / mod(delta_pos)
                tan_unit_v = norm_unit_v[::-1] * np.array((-1,1))
                u1_t = np.dot(self.vel, tan_unit_v) * tan_unit_v
                u1_n = np.dot(self.vel, norm_unit_v)
                u2_t = np.dot(p.vel, tan_unit_v) * tan_unit_v
                u2_n = np.dot(p.vel, norm_unit_v)
                v1,v2 = collision_1d(u1_n, u2_n, self.m, p.m)
                self.vel = norm_unit_v*v1 + u1_t
                p.vel    = norm_unit_v*v2 + u2_t
                self.collision_handled = 2
                p.collision_handled    = 2
    def render(self,axis):
        axis.plot(self.pos[0],self.pos[1], 'o', color = "blue")

class Membrane():
    pass

class Wall():
    pass

def run_simulation(n_water, n_solute, timesteps, temperature):
    '''
    Construct a n_timesteps x n_particles x 2 x 2 matrix, giving the position and
    velocity of each particle at each timepoint.

    Parameters
    ----------
    n_particles : TYPE
        DESCRIPTION.
    timesteps : TYPE
        DESCRIPTION.

    Returns
    -------
    np.ndarray

    '''
    tensor = np.zeros((timesteps,n_water+n_solute,2, 2))
    water = create_list_of_particles_at_temp(n_water, temperature)
    solute = create_list_of_particles_at_temp(n_solute, temperature,
                                              upper_bounds = UPPER_SOLUTE)
    particles = water + solute
    for i in range(timesteps):
        tensor[i,:,0,:] = np.array([p.pos for p in particles])
        tensor[i,:,1,:] = np.array([p.vel for p in particles])
        for p in particles: p.update(particles)
    return tensor

# def show_animation(tensor,n_water):
#     fig,ax = plt.subplots()
#     ax.set_xlim((0,1))
#     ax.set_ylim((0,1))
#     ax.vlines(0.5,0,1,linestyle="--")
#     artists = []
#     for idx,p in enumerate(tensor[0]):
#         artist, = ax.plot(p[0,0], p[0,1], 'o',
#                           color="red" if idx>n_water else "blue")
#         artists.append(artist)
#     def update(frame):
#         for artist, p in zip(artists,tensor[frame]):
#             artist.set_data(p[0,0], p[0,1])
#     global __ani
#     __ani = FuncAnimation(fig, update, interval = 50)
#     fig.show()
    
# def show_animation2(tensor,n_water,show=True):
#     fig,ax = plt.subplots()
#     ax.set_xlim((0,1))
#     ax.set_ylim((0,1))
#     ax.vlines(0.5,0,1,linestyle="--")
#     water_artist,  = ax.plot(tensor[0,:n_water,0,0], tensor[0,:n_water,0,1],'o')
#     solute_artist, = ax.plot(tensor[0,n_water:,0,0], tensor[0,n_water:,0,1],'o')
#     ratio = np.count_nonzero(tensor[0,:n_water,0,0]<0.5)/n_water
#     text = ax.text(1,1,f"Ratio of water particles = {ratio:.2f}",va="top",ha="right")
#     def update(frame):
#         water_artist.set_data(tensor[frame,:n_water,0,0], tensor[frame,:n_water,0,1])
#         solute_artist.set_data(tensor[frame,n_water:,0,0], tensor[frame,n_water:,0,1])
#         ratio = np.count_nonzero(tensor[frame,:n_water,0,0]<0.5)/n_water
#         text.set_text(f"Ratio of water particles = {ratio:.2f}")
#     global __ani
#     __ani = FuncAnimation(fig, update, interval = 5)
#     if show: fig.show()
#     return __ani

def live_animation(n_water,n_solute,temperature,show=True):
    water = create_list_of_particles_at_temp(n_water, temperature)
    solute = create_list_of_particles_at_temp(n_solute, temperature,
                                              upper_bounds = UPPER_SOLUTE)
    particles = water + solute
    fig,ax = plt.subplots()
    ax.set_xlim((0,1))
    ax.set_ylim((0,1))
    ax.vlines(0.5,0,1,linestyle="--")
    water_artist,  = ax.plot([],[],'o')
    solute_artist, = ax.plot([],[],'o')
    text = ax.text(1,1,"",va="top",ha="right")
    def update(frame):
        for p in particles: p.update(particles)
        new_pos = np.array(tuple(p.pos for p in particles))
        water_artist.set_data(new_pos[:n_water,0],new_pos[:n_water,1])
        solute_artist.set_data(new_pos[n_water:,0],new_pos[n_water:,1])
        ratio = np.count_nonzero([p.pos[0]<0.5 for p in particles[:n_water]])/n_water
        text.set_text(f"Ratio of water particles = {ratio:.2f}")
    global __ani
    __ani = FuncAnimation(fig, update, interval = 0)
    if show: fig.show()
    return __ani

# def show_animation(tensor,n_water):

#         artists.append(artist)

#     global __ani
#     __ani = FuncAnimation(fig, update, interval = 50)
#     fig.show()
    

if __name__=="__main__":
    live_animation(10,10,1)
