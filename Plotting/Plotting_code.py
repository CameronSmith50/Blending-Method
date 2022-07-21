# -*- coding: utf-8 -*-
"""
Code for plotting the blending method.
Takes in data of the form :
	data_pde_compartment_[Example Number].txt for the PDE-to-Compartment examples
	TEST[1/2/3/4]_[BLENDING/ANALYTICAL].mat for the Compartment-to-Particle examples
All parameters are defined within this script, with options for type of coupling and example number.

"""

#%% Load packages
import numpy as np
import sys
import scipy.io
import matplotlib.pyplot as plt

#%% Define the options

# Type of coupling
COUPLING_PDE_COMP = 0  # PDE-to-Compartment coupling
COUPLING_COMP_PART = 1  # Compartment-to-Particle coupling

# Check that exactly one coupling has been chosen
if COUPLING_PDE_COMP + COUPLING_COMP_PART != 1:
	sys.exit('Please choose exactly one coupling type')

# Example number
EX_1 = 1  # Pure diffusion, uniform initial condition
EX_2 = 0  # Pure diffusion, step-function initial condition
EX_3 = 0  # Morphogen gradient
EX_4 = 0  # Second-order reaction system

# Check that exactly one example has been chosen
if EX_1 + EX_2 + EX_3 + EX_4 != 1:
	sys.exit('Please choose exactly one example')

#%% Define the parameters necessary for plotting

# Parameters
x0 = 0  # Left spatial boundary
x1 = 1  # Right spatial boundary

I1 = 1/3  # Left interface position
I2 = 2/3  # Right interface position

t0 = 0  # Initial time
tf = 1  # Final time

dt = 1e-2  # Recording time-step
nt = tf/dt + 1  # Number of recording steps

hc = 1/30  # Compartment width
K = 30  # Number of compartments

hp = 1/300  # PDE mesh width

# Temporal mesh
t_MESH = np.arange(t0,tf+dt,dt)

# Spatial mesh for PDE-to-Compartment
if COUPLING_PDE_COMP == 1:
	LEFT_MESH = np.arange(x0,I1,hp)
	BLEND_MESH = np.arange(I1+hc/2,I2,hc)
	RIGHT_MESH = np.arange(I2+hc/2,x1,hc)
	AN_MESH = np.arange(x0,x1+hp,hp)

# Spatial mesh for Compartment-to-Particle
else:
	LEFT_MESH = np.arange(x0+hc/2,I1,hc)
	BLEND_MESH = np.arange(I1+hc/2,I2,hc)
	RIGHT_MESH = np.arange(I2+hc/2,x1,hc)
	AN_MESH = np.arange(x0,x1+hp,hp)


#%% Import data and create strings for saving

# PDE-to-Compartment
if COUPLING_PDE_COMP == 1:
	COUPLE_STR = 'pde_comp'
	if EX_1 == 1:
		BLENDING = np.load('./Data_hybrid_compartment_pde/experiment/TEST_PROBLEM_1.npy')
		ANALYTICAL = scipy.io.loadmat('./Data_hybrid_compartment_particle/TEST_PROBLEM_1/TEST1_ANALYTICAL')['u']
		EX_STR = '1'
	elif EX_2 == 1:
		BLENDING = np.load('./Data_hybrid_compartment_pde/experiment/TEST_PROBLEM_2.npy')
		ANALYTICAL = scipy.io.loadmat('./Data_hybrid_compartment_particle/TEST_PROBLEM_2/TEST2_ANALYTICAL')['u']
		EX_STR = '2'
	elif EX_3 == 1:
		BLENDING = np.load('./Data_hybrid_compartment_pde/experiment/TEST_PROBLEM_3.npy')
		ANALYTICAL = scipy.io.loadmat('./Data_hybrid_compartment_particle/TEST_PROBLEM_3/TEST3_ANALYTICAL')['u']
		EX_STR = '3'
	else:
		BLENDING = np.load('./Data_hybrid_compartment_pde/experiment/TEST_PROBLEM_4.npy')
		ANALYTICAL = np.load('./Data_hybrid_compartment_pde/experiment/TEST_4_GROUND_TRUTH_500SIM.npy')
		EX_STR = '4'

		# Update parameter values
		x1 = 10
		I1 = 10/3
		I2 = 20/3
		hc = 1/3
		hp = 1/30
		LEFT_MESH = np.arange(x0+hp/2,I1,hp)
		BLEND_MESH = np.arange(I1+hc/2,I2,hc)
		RIGHT_MESH = np.arange(I2+hc/2,x1,hc)
		AN_MESH = np.arange(x0+hc/2,x1,hc)

else:
	COUPLE_STR = 'comp_part'
	if EX_1 == 1:
		BLENDING = scipy.io.loadmat('./Data_hybrid_compartment_particle/TEST_PROBLEM_1/TEST1_BLENDING')['meso_micro']
		ANALYTICAL = scipy.io.loadmat('./Data_hybrid_compartment_particle/TEST_PROBLEM_1/TEST1_ANALYTICAL')['u']
		EX_STR = '1'
	elif EX_2 == 1:
		BLENDING = scipy.io.loadmat('./Data_hybrid_compartment_particle/TEST_PROBLEM_2/TEST2_BLENDING')['meso_micro']
		ANALYTICAL = scipy.io.loadmat('./Data_hybrid_compartment_particle/TEST_PROBLEM_2/TEST2_ANALYTICAL')['u']
		EX_STR = '2'
	elif EX_3 == 1:
		BLENDING = scipy.io.loadmat('./Data_hybrid_compartment_particle/TEST_PROBLEM_3/TEST3_BLENDING')['MG']
		ANALYTICAL = scipy.io.loadmat('./Data_hybrid_compartment_particle/TEST_PROBLEM_3/TEST3_ANALYTICAL')['u']
		EX_STR = '3'
	else:
		BLENDING = scipy.io.loadmat('./Data_hybrid_compartment_particle/TEST_PROBLEM_4/TEST4_blend')['second_order_blend']
		ANALYTICAL = scipy.io.loadmat('./Data_hybrid_compartment_particle/TEST_PROBLEM_4/TEST4_particle')['second_order_particle']
		EX_STR = '4'

		# Update parameter values
		x1 = 10
		I1 = 10/3
		I2 = 20/3
		hc = 1/3
		LEFT_MESH = np.arange(x0+hc/2,I1,hc)
		BLEND_MESH = np.arange(I1+hc/2,I2,hc)
		RIGHT_MESH = np.arange(I2+hc/2,x1,hc)
		AN_MESH = np.arange(x0+hc/2,x1,hc)

#%% Calculate relative errors

# The number of PDE-Compartment time-steps to Analytical time-steps
if EX_1 == 1 or EX_2 == 1 or EX_4 == 1:
	factor = 100
elif EX_3 == 1 :
	factor = 100

# PDE-to-Compartment
if COUPLING_PDE_COMP == 1:
	if EX_4 == 1:
		ANALYTICAL_PARTS_LEFT = np.sum(ANALYTICAL[1:101,0:10],axis=1)
	else:
		ANALYTICAL_PARTS_LEFT = (np.sum(ANALYTICAL[1:101,0:100],1)+np.sum(ANALYTICAL[1:101,1:101],1))*hp/2
	REL_ERROR_LEFT = (np.sum(BLENDING[np.arange(factor,factor*(nt-1)+1,factor,dtype=int),0:100],axis=1)*hp - ANALYTICAL_PARTS_LEFT)/ANALYTICAL_PARTS_LEFT

	if EX_4 == 1:
		ANALYTICAL_PARTS_BLEND = np.sum(ANALYTICAL[1:101,10:20],axis=1)
	else:
		ANALYTICAL_PARTS_BLEND = (np.sum(ANALYTICAL[1:101,100:200],1)+np.sum(ANALYTICAL[1:101,101:201],1))*hp/2
	REL_ERROR_BLEND = (np.sum(BLENDING[np.arange(factor,factor*(nt-1)+1,factor,dtype=int),200:210],axis=1)*hc - ANALYTICAL_PARTS_BLEND)/ANALYTICAL_PARTS_BLEND

	if EX_4 == 1:
		ANALYTICAL_PARTS_RIGHT = np.sum(ANALYTICAL[1:101,20:30],axis=1)
	else:
		ANALYTICAL_PARTS_RIGHT = (np.sum(ANALYTICAL[1:101,200:300],1)+np.sum(ANALYTICAL[1:101,201:301],1))*hp/2
	REL_ERROR_RIGHT = (np.sum(BLENDING[np.arange(factor,factor*(nt-1)+1,factor,dtype=int),210:220],axis=1)*hc - ANALYTICAL_PARTS_RIGHT)/ANALYTICAL_PARTS_RIGHT

# Compartment-to-particle
else:
	if EX_4 == 1:
		ANALYTICAL_PARTS_LEFT = np.sum(ANALYTICAL[1:101,0:10],axis=1)
	else:
		ANALYTICAL_PARTS_LEFT = 0.5*(np.sum(ANALYTICAL[1:101,0:100],1)+np.sum(ANALYTICAL[1:101,1:101],1))*hp
	REL_ERROR_LEFT = (np.sum(BLENDING[1:101,0:10],1) - ANALYTICAL_PARTS_LEFT)/ANALYTICAL_PARTS_LEFT

	if EX_4 == 1:
		ANALYTICAL_PARTS_BLEND = np.sum(ANALYTICAL[1:101,10:20],axis=1)
	else:
		ANALYTICAL_PARTS_BLEND = 0.5*(np.sum(ANALYTICAL[1:101,100:200],1)+np.sum(ANALYTICAL[1:101,101:201],1))*hp
	REL_ERROR_BLEND = (np.sum(BLENDING[1:101,10:20],1) - ANALYTICAL_PARTS_BLEND)/ANALYTICAL_PARTS_BLEND

	if EX_4 == 1:
		ANALYTICAL_PARTS_RIGHT = np.sum(ANALYTICAL[1:101,20:30],axis=1)
	else:
		ANALYTICAL_PARTS_RIGHT = 0.5*(np.sum(ANALYTICAL[1:101,200:300],1)+np.sum(ANALYTICAL[1:101,201:301],1))*hp
	REL_ERROR_RIGHT = (np.sum(BLENDING[1:101,20:30],1) - ANALYTICAL_PARTS_RIGHT)/ANALYTICAL_PARTS_RIGHT

# Find the y limits
y_MIN = np.min((np.min(REL_ERROR_LEFT[1:101:2]),np.min(REL_ERROR_BLEND[1:101:2]),np.min(REL_ERROR_RIGHT[1:101:2])))
y_MAX = np.max((np.max(REL_ERROR_LEFT[1:101:2]),np.max(REL_ERROR_BLEND[1:101:2]),np.max(REL_ERROR_RIGHT[1:101:2])))

MAXIMUM = np.max((np.abs(y_MIN),np.abs(y_MAX)))

# Give the limits
y0 = -MAXIMUM
y1 = MAXIMUM

# Set fontsize
fsize = 100 #XXX
ticksize = 50

#%% Plotting

# Times to plot at
plot_times = np.array([0,tf/10,tf])
if EX_3 == 1:
	plot_times[1] = tf/100
plot_inds_coarse = np.array(plot_times/dt,dtype=int)

# Set the y limits for the different examples
if EX_1 == 1:
	yLimUpperDens = 1200
	yLimUpperErr = 0.008
elif EX_2 == 1:
	yLimUpperDens = 1800
	yLimUpperErr = 0.01
elif EX_3 == 1:
	yLimUpperDens = 3500
	yLimUpperErr = 0.015
elif EX_4 == 1:
	yLimUpperDens = 60
	yLimUpperErr = 0.03

# Index for plotting references
j = 1

# Loop through the plot times
for i in plot_inds_coarse:

	# Start a figure
	plt.figure(j, figsize=(24,12.5))  #XXX Start figure with indentifier i
	plt.xlabel('x', fontsize=fsize)  # xlabel
	plt.ylabel('Particle density', fontsize=fsize)  # ylabel
#	plt.title('Time = %g' % t_MESH[i], fontsize=fsize)
	plt.rc('xtick', labelsize=fsize)
	plt.rc('ytick', labelsize=fsize)
	plt.xlim((x0,x1))


	# Split into two cases depending on the coupling
	# PDE-to-Compartment
	if COUPLING_PDE_COMP == 1:

		# Calculate maximum y value
		max_y = max(np.max(BLENDING),np.max(ANALYTICAL))*1.2

		plt.plot(LEFT_MESH,BLENDING[factor*i,0:100], lw=4, color='#009900')  # PDE
		plt.bar(BLEND_MESH, BLENDING[factor*i,200:210], 7*hc/8, color='#FF0000', edgecolor='black')  # Blended or #004C80
		plt.bar(RIGHT_MESH, BLENDING[factor*i,210:220], 7*hc/8, color='#0000FF', edgecolor='black')  # Compartments
		if EX_4 == 1:
			plt.plot(AN_MESH, ANALYTICAL[i,]/hc, ls='--', lw=4, color='#000000')  # Analytical
		else:
			plt.plot(AN_MESH, ANALYTICAL[i,], ls='--', lw=4, color='#000000')  # Analytical

	# Compartment-to-Particle
	else:

		# Calculate maximum y value
		max_y = max(np.max(BLENDING/hc),np.max(ANALYTICAL))*1.2

		plt.bar(LEFT_MESH, BLENDING[i,0:10]/hc, 7*hc/8, color='#0000FF', edgecolor='black')  # Compartments
		plt.bar(BLEND_MESH, BLENDING[i,10:20]/hc, 7*hc/8, color='#FF0000', edgecolor='black')  # Blended or #808080
		plt.bar(RIGHT_MESH, BLENDING[i,20:30]/hc, 7*hc/8, color='#FFFF00', edgecolor='black')  # Particles
		if EX_4 == 1:
			plt.plot(AN_MESH, ANALYTICAL[i,]/hc, ls='--', lw=4, color='#000000')  # Analytical
		else:
			plt.plot(AN_MESH, ANALYTICAL[i,], ls='--', lw=4, color='#000000')  # Analytical

	# Add in the common parts
	plt.ylim((0,max_y))
	ax = plt.gca()
	plt.ylim(0,yLimUpperDens)
	plt.plot(I1*np.ones((2,1)), np.array([0,yLimUpperDens]), lw=4, color='#FF0000')  # Interface 1
	plt.plot(I2*np.ones((2,1)), np.array([0,yLimUpperDens]), lw=4, color='#FF0000')  # Interface 2
	plt.xticks(fontsize=ticksize)
	plt.yticks(fontsize=ticksize)

	# Move figure in frame
	pos1 = ax.get_position()
	pos2 = [pos1.x0 + 0.05, pos1.y0 + 0.06,  pos1.width, pos1.height]
	ax.set_position(pos2)
	ax.tick_params(axis='x', pad=15)

	# Save the image
	plt.savefig('./Figures/coupling_%s_example_%s_time_%d.eps' % (COUPLE_STR,EX_STR,int(100*t_MESH[i])))

	# Update the dummy index
	j += 1

#%% Errors

# # Left-hand side
# plt.figure(4, figsize=(25,14.5))  #XXX Start figure with indentifier i+1
# plt.xlabel('Time', fontsize=fsize)  # xlabel
# plt.ylabel('Relative error', fontsize=fsize)  # ylabel
# plt.xlim((t0,tf))
# plt.ylim((-yLimUpperErr,yLimUpperErr))
# #if COUPLING_PDE_COMP == 1:
# #	plt.title('PDE Subdomain', fontsize=fsize)
# #else:
# #	plt.title('Compartment Subdomain', fontsize=fsize)

# # Plot 0 in a red dashed line
# plt.plot(np.array([t0,tf]),np.zeros((2,1)), lw=2, ls='--', color='#FF0000')

# # Plot the error in black
# plt.plot(t_MESH[1:101], REL_ERROR_LEFT, lw=3, color='#000000')

# # Save the figure
# plt.savefig('./Figures/coupling_%s_example_%s_rel_error_left.eps' % (COUPLE_STR,EX_STR), format='eps', dpi=1000)

# # Blended region
# plt.figure(5, figsize=(25,14.5))  #XXX Start figure with indentifier i+1
# plt.xlabel('Time', fontsize=fsize)  # xlabel
# plt.ylabel('Relative error', fontsize=fsize)  # ylabel
# plt.xlim((t0,tf))
# plt.ylim((-yLimUpperErr,yLimUpperErr))
# #plt.title('Blended Subdomain', fontsize=fsize)

# # Plot 0 in a red dashed line
# plt.plot(np.array([t0,tf]),np.zeros((2,1)), lw=2, ls='--', color='#FF0000')

# # Plot the error in black
# plt.plot(t_MESH[1:101], REL_ERROR_BLEND, lw=3, color='#000000')

# # Save the figure
# plt.savefig('./Figures/coupling_%s_example_%s_rel_error_blend.eps' % (COUPLE_STR,EX_STR), format='eps', dpi=1000)

# # Right-hand side
# plt.figure(6, figsize=(25,14.5))  #XXX Start figure with indentifier i+1
# plt.xlabel('Time', fontsize=fsize)  # xlabel
# plt.ylabel('Relative error', fontsize=fsize)  # ylabel
# plt.xlim((t0,tf))
# plt.ylim((-yLimUpperErr,yLimUpperErr))
# #if COUPLING_PDE_COMP == 1:
# #	plt.title('Compartment Subdomain', fontsize=fsize)
# #else:
# #	plt.title('Particle Subdomain', fontsize=fsize)

# # Plot 0 in a red dashed line
# plt.plot(np.array([t0,tf]),np.zeros((2,1)), lw=2, ls='--', color='#FF0000')

# # Plot the error in black
# plt.plot(t_MESH[1:101], REL_ERROR_RIGHT, lw=3, color='#000000')

# # Save the figure
# plt.savefig('./Figures/coupling_%s_example_%s_rel_error_right.eps' % (COUPLE_STR,EX_STR), format='eps', dpi=1000)

# All three in one
# Right-hand side
ax = plt.figure(7, figsize=(24,12.5))  #XXX Start figure with indentifier i+1
plt.xlabel('Time', fontsize=fsize)  # xlabel
plt.ylabel('Relative error', fontsize=fsize)  # ylabel
plt.xlim((t0,tf))
plt.ylim((-yLimUpperErr,yLimUpperErr))
#plt.title('Relative errors', fontsize=fsize)

# Plot 0 in a red dashed line
plt.plot(np.array([t0,tf]),np.zeros((2,1)), lw=2, ls='--', color='#000000')

# Plot the error in black
if COUPLING_PDE_COMP == 1:
 	left, = plt.plot(t_MESH[1:101:2], REL_ERROR_LEFT[0:101:2], '--', lw=2, color='#009900', label='PDE subdomain')
 	blend, = plt.plot(t_MESH[1:101:2], REL_ERROR_BLEND[0:101:2], '-', lw=2, color='#FF0000', label='Blended subdomain')
 	right, = plt.plot(t_MESH[1:101:2], REL_ERROR_RIGHT[0:101:2], ':', lw=2, color='#0000FF', label='Compartment subdomain')
else:
 	left, = plt.plot(t_MESH[1:101:2], REL_ERROR_LEFT[0:101:2], '--', lw=2, color='#0000FF', label='Compartment subdomain')
 	blend, = plt.plot(t_MESH[1:101:2], REL_ERROR_BLEND[0:101:2], '-', lw=2, color='#FF0000', label='Blended subdomain')
 	right, = plt.plot(t_MESH[1:101:2], REL_ERROR_RIGHT[0:101:2], ':', lw=2, color='#FF9900', label='Particle subdomain')

if EX_1 == 1 or EX_3 == 1:
	plt.legend(loc='lower left', fontsize = 48, borderaxespad=0, borderpad = 0.1, fancybox=False)
elif EX_2 == 1:
	plt.legend(loc='lower right', fontsize = 48, borderaxespad=0, borderpad = 0.1, fancybox=False)
else:
	plt.legend(loc='upper right', fontsize = 48, borderaxespad=0, borderpad = 0.1, fancybox=False)
plt.xticks(fontsize=ticksize)
plt.yticks(fontsize=ticksize)
ax = plt.gca()
pos1 = ax.get_position()
pos2 = [pos1.x0 + 0.06, pos1.y0 + 0.07,  pos1.width, pos1.height]
ax.set_position(pos2)
ax.tick_params(axis='x', pad=15)

plt.savefig('./Figures/coupling_%s_example_%s_rel_error_all.eps' % (COUPLE_STR, EX_STR), format='eps', dpi=1000)
