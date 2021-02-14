# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 16:29:44 2021

@author: Vivian Imbriotis
"""

import osmosisbackend as osb
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def live_animation():
    fig,ax = plt.subplots()
    nwater = 100
    state = osb.construct_state(200,nwater,1)
    water_artist, = ax.plot([],[],'o')
    solute_artist,= ax.plot([],[],'o')
    def update(frame):
        nonlocal state
        positions = osb.get_pos_from_state(state)
        water_artist.set_data(positions[::2][:nwater],positions[1::2][:nwater])
        solute_artist.set_data(positions[::2][nwater:],positions[1::2][nwater:])
        state = osb.update_state(state)
    global __ani
    __ani = FuncAnimation(fig,update, interval=0)
    fig.show()
    return __ani

if __name__=="__main__":
    ani = live_animation()
