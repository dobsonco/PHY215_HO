import numpy as np
import matplotlib.pyplot as plt
import random as rd
from collections import Counter
import streamlit as st

class Nernst():

    def find_potential(self,T,z,Xi,Xo):
        """
        T: Temperature in Kelvin
        z: Valence electrons of ion
        Xi: Inner membrane concentration (mM)
        Xo: Outer membrane concentration (mM)
        returns: Membrane potential for single species
        """
        membranePotential = (8.314*T)/(z*96485)*np.log(Xo/Xi)
        return membranePotential
    
    def plot_potential(self,volts,name,color):
        '''
        This function plots the potential across a membrane
        volts: Float, The voltage potential calculated by find_potential
        name: String, The name of the ionic species
        '''
        points = [0,1,2,3]
        values= [0,0,volts,volts]
        plt.scatter(points,values,label=name,color=color)
        plt.plot(points,values,color=color)
        plt.ylim(-100,100)
        plt.xlim(0,3)
        xlabs = np.arange(0, 4, 1.0)
        x_tick_names = ['','Channel Closed', 'Channel Open','']
        plt.xticks(xlabs,x_tick_names,color='lightgray')
        plt.text(2, volts+3, str(round(volts, 2)))
        plt.title(f'Membrane Potential for {name}',color='lightgray')
        plt.xlabel("Membrane Permeability",color='lightgray')
        plt.ylabel("Voltage (mV)",color='lightgray')
        plt.legend()
        
    def plot_positions(self,ion_object_list):
        '''
        This function takes in a list of ion objects, and plots them on a
        scatterplot. All ions of the same name are given the same color, but
        said color is random. Also returns a variable (called unique_colors)
        with a list of colors corresponding to the different ionic species in solution.
        '''
        names = []
        colors = []
        [names.append(x.get_name()) for x in ion_object_list]
        counted = dict(Counter(names))
        for key in counted:
            color = "#"+''.join([rd.choice('0123456789ABCDEF') for j in range(6)])
            for ion in range(counted[key]):
                colors.append(color)

        unique_colors = []
        [unique_colors.append(x) for x in colors if x not in unique_colors]
        
        for i,obj in enumerate(ion_object_list):
            position = obj.get_position()
            if obj.get_name() == ion_object_list[i-1].get_name():
                plt.scatter(position[0],position[1],color=colors[i])
            else:
                plt.scatter(position[0],position[1],color=colors[i],label=names[i])
        plt.plot([0,0],[0,1.1],c='salmon',linewidth=2,label='Membrane')
        
        xlabs = np.arange(-1.0,1.0,0.5)
        x_tick_names = ['','In','','Out']
        plt.xticks(xlabs,x_tick_names,color='lightgray')
        
        plt.xlim(-1,1)
        plt.ylim(0,1.01)
        plt.title('Distrubution of Ions',color='lightgray')
        plt.legend()
        plt.tick_params(left=False,bottom=True,labelleft=False,labelbottom=True)
        return unique_colors
        
    def setup(self,mol_dict,total_particles=100):
        '''
        This function takes in a dictionary with ionic species as the key,
        the valence electrons in the first column, concentration inside the cell
        in the second, and the concentration outside the cell in the third and returns
        a list of objects that is proportional to the concentrations. The returned list
        also has all of the correct attributes such that they can be easily plotted.
        Has no effect on the potential.
        '''
        total_particles = total_particles
        ion_object_list = []

        total = 0
        for key in mol_dict:
            for value in mol_dict[key][1:]:
                total += value

        chunk_size = total_particles / total
        
        for key in mol_dict:
            for i,dont_use in enumerate(mol_dict[key][1:]):
                num_ions = int(round(mol_dict[key][i+1]*chunk_size,0))
                for j in range(num_ions):
                    if i%2 != 0:
                        new_ion = ion(key,mol_dict[key][0],'out')
                        ion_object_list.append(new_ion)
                    else:
                        new_ion = ion(key,mol_dict[key][0],'in')
                        ion_object_list.append(new_ion)
        return(ion_object_list)
    
Nernst = Nernst()

class ion(object):
    
    def __init__(self,ion_name,valence,IO):
        """ Initialize ion object.
        Inputs are
        ion_name: the type of ion (string)
        valence: number of valence electrons (int)
        location: string, either in or out
        """
        self.ion_name = ion_name
        self.valence = valence
        self.IO = IO
        if IO == 'in':
            self.position = (rd.uniform(-0.99,-0.02),rd.uniform(0.01,0.99))
        else:
            self.position = (rd.uniform(0.02,0.99),rd.uniform(0.01,0.99))
        
    def get_name(self):
        return(self.ion_name)
    
    def get_valence(self):
        return(self.valence)
    
    def get_IO(self):
        return(self.IO)
    
    def get_position(self):
        return(self.position)

if __name__ == '__main__':
    
    molarity_dict = {
#                  z  in  out
        'sodium': (1, 0.1, 2),
        'calcium': (2, 4, 0.1),
        'potassium': (1, 1, 4),
    }
    
    ion_object_list = Nernst.setup(molarity_dict)
    
##################################################### Plot Ions
    figure = plt.figure(figsize=(10,10))
    plt.subplot(2,1,2)
    unique_colors = Nernst.plot_positions(ion_object_list)
    
##################################################### Plot Potential

#Add code here that takes in a dictionary of colors for the different species
#so that the colors in plot one and two match
    plt.subplot(2,1,1)
    temp = 300
    
    for i,key in enumerate(molarity_dict):
        z = molarity_dict[key][0]
        inside = molarity_dict[key][1]
        outside = molarity_dict[key][2]
        potential = Nernst.find_potential(temp,z,inside,outside)
        Nernst.plot_potential(potential*1000,key,unique_colors[i])

    plt.style.use("dark_background")
    plt.tight_layout()
    
    st.pyplot(figure)