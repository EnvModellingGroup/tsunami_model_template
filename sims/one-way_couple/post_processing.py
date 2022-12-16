import shutil
from thetis import *
import os.path
import sys
import math
import sys
import os.path
import numpy as np
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import params


# where would you like to store the output of this analysis?
output_dir = 'analysis'

# create BSS per timestep? Probbaly false unless you're doing tsunamis...
BSS_TS = True

# You *MAY* need to edit below this line
# Make sure below matches your main run file as much as possible
# *if* anything goes wrong with the analysis
#============================================================#
mesh2d = Mesh(os.path.join(os.path.pardir,os.path.pardir,params.mesh_file))

# How long does your simulations run for (s)
t_end = params.end_time # which is the start file?
start_file = 0
# how often are exports produced in the main run?
t_export = params.output_time

# where are your thetis output files (note, do not include the hdf5 directory)
thetis_dir = params.output_dir

t_n = int((t_end/t_export) - start_file + 1)
thetis_times = t_export*np.arange(t_n) + t_export

P1 = FunctionSpace(mesh2d, "CG", 1)
# read bathymetry code
chk = DumbCheckpoint('bathymetry', mode=FILE_READ)
bathymetry2d = Function(P1)
chk.load(bathymetry2d,  name='bathymetry')
chk.close()
chk = DumbCheckpoint('manning', mode=FILE_READ)
manning = Function(bathymetry2d.function_space(), name='manning')
chk.load(manning)
chk.close()

# --- create solver ---
solverObj = solver2d.FlowSolver2d(mesh2d, bathymetry2d)
options = solverObj.options
options.use_nonlinear_equations = True
options.simulation_export_time = t_export
options.simulation_end_time = t_end
options.output_directory =  thetis_dir
options.check_volume_conservation_2d = True
options.fields_to_export = ['uv_2d', 'elev_2d']
options.manning_drag_coefficient = Constant(1.0)
options.horizontal_viscosity = Constant(1.0)
options.coriolis_frequency = Constant(1.0)
options.timestep = 1.0
options.use_grad_div_viscosity_term = True
options.use_wetting_and_drying = True
options.wetting_and_drying_alpha = Constant(1.0)
options.element_family = "dg-dg"
options.swe_timestepper_type = 'CrankNicolson'
options.swe_timestepper_options.implicitness_theta = 1.0
options.swe_timestepper_options.use_semi_implicit_linearization = True
options.use_grad_div_viscosity_term = True
options.use_grad_depth_viscosity_term = False

# set boundary/initial conditions code
solverObj.bnd_functions['shallow_water'] = {
    666: {'elev': 0.0},
    1000: {'un': 0.0},
    #set closed boundaries to zero velocity
}

# we need bathy and manning on the same mesh as the elev and vel
P1DG = FunctionSpace(mesh2d, "DG", 1)
bathydg = project(bathymetry2d, P1DG)
manningdg = project(manning, P1DG)


uv = Function(P1DG, name='vel_2d')
u_data_set = np.empty((t_n, uv.dat.data.shape[0]))
v_data_set = np.empty((t_n, uv.dat.data.shape[0]))
elev_data_set = np.empty((t_n, uv.dat.data.shape[0])) # elevation is also DG, so we can just do this.
count = 0
if BSS_TS:
    bss = Function(P1DG, name="BSS")
    bathy = bathydg.dat.data[:]
    bss_file = File( output_dir + '/bss.pvd')
for i in range(start_file,int(t_end/t_export)+1):
    solverObj.load_state(i)
    u_data_set[count, :] = solverObj.fields.uv_2d.dat.data[:,0]
    v_data_set[count, :] = solverObj.fields.uv_2d.dat.data[:,1]
    elev_data_set[count, :] = solverObj.fields.elev_2d.dat.data[:]
    if BSS_TS:
        speed = np.sqrt(u_data_set[count, :]*u_data_set[count, :] + v_data_set[count, :]*v_data_set[count, :])
        elev_bathy = elev_data_set[count, :] + bathy
        man = manningdg.dat.data[:]
        tau_b = np.array(1024*9.81*man*man*speed*speed / (elev_bathy)**(1./3.))
        tau_b[ elev_bathy < 0.001] = 0.0 # we have < 1mm of water
        tau_b[ tau_b < 0.0 ] = 0.0 # we had no water (shouldn't happen due to above, but just in case)
        bss.dat.data[:] = tau_b 
        bss_file.write(bss)
    count += 1

ave_speed = [] # average over speeds
max_speed = [] # max over speeds
ave_bss = [] # ave of bss calc
max_bss = [] # max of bss calc
ave_vel = [] # vector of ave u and ave v
max_vel = [] # vector of when max speed occurs

for i in range(uv.dat.data.shape[0]): # loop over nodes in the Function mesh
    u_vel = np.array(u_data_set[:,i]) # timeseries of u, v and elevation
    v_vel = np.array(v_data_set[:,i])
    elev = np.array(elev_data_set[:,i])
    speed = np.sqrt(u_vel*u_vel + v_vel*v_vel)
    ave_speed.append(np.mean(speed))
    max_speed.append(np.max(speed))

    man = np.array(manningdg.dat.data[i])
    bathy = np.array(bathydg.dat.data[i])
    elev_bathy = elev+bathy
    power = np.sign(elev_bathy) * (np.abs(elev_bathy)) ** (1. / 3.)
    tau_b = np.array(1024*9.81*man*man*speed*speed / power)
    tau_b[ elev_bathy < 0.01] = 0.0
    tau_b[ tau_b < 0.0 ] = 0.0
    ave_bss.append(np.mean(tau_b))
    max_bss.append(np.max(tau_b))

    ave_vel.append([np.mean(u_vel), np.mean(v_vel)])
    max_vel.append([u_vel[np.argmax(speed)],v_vel[np.argmax(speed)]])


avespeed = Function(P1DG, name="AveSpeed")
avespeed.dat.data[:] = np.array(ave_speed)
File( output_dir + '/average_speed.pvd').write(avespeed)
chk = DumbCheckpoint(output_dir + '/average_speed', mode=FILE_CREATE)
chk.store(avespeed, name='AveSpeed')
maxspeed = Function(P1DG, name="MaxSpeed")
maxspeed.dat.data[:] = np.array(max_speed)
File( output_dir + '/max_speed.pvd').write(maxspeed)
chk = DumbCheckpoint(output_dir + '/max_speed', mode=FILE_CREATE)
chk.store(maxspeed, name='MaxSpeed')
avebss = Function(P1DG, name="AveBSS")
avebss.dat.data[:] = np.array(ave_bss)
File( output_dir + '/average_bss.pvd').write(avebss)
chk = DumbCheckpoint(output_dir + '/average_bss', mode=FILE_CREATE)
chk.store(avebss, name='AveBSS')
maxbss = Function(P1DG, name="MaxBSS")
maxbss.dat.data[:] = np.array(max_bss)
File( output_dir + '/max_bss.pvd').write(maxbss)
chk = DumbCheckpoint(output_dir + '/max_bss', mode=FILE_CREATE)
chk.store(maxbss, name='MaxBSS')
# now the vectors
P1DG = VectorFunctionSpace(mesh2d, "DG", 1)
avevel = Function(P1DG, name='ave_vel')
avevel.dat.data[:] = np.array(ave_vel)
File( output_dir + '/average_vel.pvd').write(avevel)
chk = DumbCheckpoint(output_dir + '/average_vel', mode=FILE_CREATE)
chk.store(avevel, name='AveVel')
maxvel = Function(P1DG, name='max_vel')
maxvel.dat.data[:] = np.array(max_vel)
File( output_dir + '/max_vel.pvd').write(maxvel)
chk = DumbCheckpoint(output_dir + '/max_vel', mode=FILE_CREATE)
chk.store(maxvel, name='MaxVel')


