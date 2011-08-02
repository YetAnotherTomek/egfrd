#!/usr/bin/env python
# -*- coding: utf-8 -*-


LOGGER=True
VTK_LOGGER=True


import sys
from egfrd import *
from bd import *
import model
import gfrdbase
import myrandom


def singlerun(T):


    ### Define system constants
    L=1e-3
    matrix_size = 3
    N_particles = 24
    
    sigma = 3e-5
    r0 = sigma
    D = 2e-12
    D_tot = D
    v = 1
    kf = 10 * sigma * D_tot

    radius_cyl = 5e-5


    ### Create the model and geometric structures
    m = model.ParticleModel(L)
    
    cyl = model.create_cylindrical_surface('cyl', [L/2, 0, L/2], radius_cyl, [0, 1, 0], L*1.0)
    m.add_structure(cyl)

    pla = model.create_planar_surface('pla', [0,0,0], [0,1,0], [1,0,0], L, L)
    m.add_structure(pla)

    A = model.Species('A', D, sigma/2)
    m.add_species_type(A)
    Acyl = model.Species('Acyl', D/10, sigma/2, 'cyl', v)
    m.add_species_type(Acyl)
    Apla = model.Species('Apla', D/10, sigma/2, 'pla')
    m.add_species_type(Apla)
    
    B = model.Species('B', D, sigma/2)
    m.add_species_type(B)
    C = model.Species('C', D, sigma/2)
    m.add_species_type(C)

    ruleABC = model.create_binding_reaction_rule(A, B, C, kf)
    m.network_rules.add_reaction_rule(ruleABC)
    
    print "Model and geometric structures defined."


    ### Create the world
    w = gfrdbase.create_world(m, matrix_size)
    nrw = gfrdbase.create_network_rules_wrapper(m)
    s = EGFRDSimulator(w)
    s.bd_dt_factor = 1

    ### Place particles    
    gfrdbase.throw_in_particles(w,A,20)	# typically 20
    gfrdbase.throw_in_particles(w,Acyl,5)
    gfrdbase.throw_in_particles(w,Apla,5)
    
    print "World created and particles placed."
    
    
    ### Define and initialize the loggers
    if VTK_LOGGER == True:
        from visualization import vtklogger
        # Write vtk files. See vtklogger.py. Use VTK or Paraview to 
        # visualize.  
        vtklogger = vtklogger.VTKLogger(s, 'data/run', extra_particle_step=True)

    if LOGGER == True:
        from logger import Logger
        # Write log files. See logger.py. Use visualizer.py to visualize.
        l = Logger('example')
        l.start(s)
        #l.set_particle_out_interval(1e-4) #1e-3)
        
    print "Loggers initialized."


    ### RUN IT
    print "Starting the simulation."
    end_time = T
    s.step()

    while 1:
        next_time = s.get_next_time()
        print next_time
       	vtklogger.log()
        
        if next_time > end_time:
            s.stop(end_time)
            if VTK_LOGGER == True:
                vtklogger.stop()
                
            print "Simulation end."
            break
        
        s.step()
        if s.last_reaction:
            print 'reaction'
            return 0.0, s.t
    

''' MAIN

'''
if __name__ == '__main__':

    singlerun(500.0)
