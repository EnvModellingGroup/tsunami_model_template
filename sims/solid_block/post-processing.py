import shutil
from thetis import *
import os.path
import sys
import math
import numpy as np
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import params

# where should the output of this analysis go
output_dir = 'analysis'
create_directory(output_dir)

# where is the output of your model?
thetis_dir = params.output_dir

# was this run created with the DumbCheckpoint code? If so, make this True
legacy_run = False

# You *MAY* need to edit below this line
# Make sure below matches your main run file as much as possible
# *if* anything goes wrong with the analysis
#============================================================#

# making an assumption here on what the hdf5 output is called
chk = CheckpointFile("output/hdf5/Elevation2d_00000.h5",'r')
thetis_mesh = chk.load_mesh()

chk = CheckpointFile('bathymetry.h5','r')
bathymetry2d = chk.load_function(thetis_mesh,'bathymetry')
chk.close()
chk = CheckpointFile('manning.h5','r')
manning = chk.load_function(thetis_mesh, 'manning')
chk.close()

# How long does your simulations run for (s)
t_end = params.end_time #40 days (i.e. 30 days of analysis)
# how often are exports produced in the main run?
t_export = params.output_time


t_n = int((t_end - t_start) / t_export) + 1
thetis_times = t_start + t_export*np.arange(t_n)


P1 = FunctionSpace(thetis_mesh, "CG", 1)

# --- create solver ---
solverObj = solver2d.FlowSolver2d(thetis_mesh, Constant(10))
options = solverObj.options
options.simulation_export_time = t_export
options.simulation_end_time = t_end
options.output_directory =  thetis_dir
options.manning_drag_coefficient = Constant(1.0)
options.horizontal_viscosity = Constant(1.0)
options.element_family = "dg-dg"

# we need bathy and manning on the same mesh as the elev and vel
P1DG = FunctionSpace(thetis_mesh, "DG", 1)
manningdg = project(manning, P1DG)
bathydg = project(bathymetry2d, P1DG)

uv = Function(P1DG, name='vel_2d')
u_data_set = np.empty((t_n, uv.dat.data.shape[0]))
v_data_set = np.empty((t_n, uv.dat.data.shape[0]))
elev = Function(P1DG, name='elev_2d')
elev_data_set = np.empty((t_n, elev.dat.data.shape[0]))
count = 0
bss = Function(P1DG, name="BSS")
bathy = bathydg.dat.data[:]
bss_file = File( output_dir + '/bss.pvd')

for t in thetis_times:
    PETSc.Sys.Print('Reading h5 files. Time ',int(t/t_export),t)
    solverObj.load_state(int(t/t_export), legacy_mode=legacy_run)
    u_data_set[count, :] = solverObj.fields.uv_2d.dat.data[:,0]
    v_data_set[count, :] = solverObj.fields.uv_2d.dat.data[:,1]
    elev_data_set[count, :] = solverObj.fields.elev_2d.dat.data[:]
    speed = np.sqrt(u_data_set[count, :]*u_data_set[count, :] + v_data_set[count, :]*v_data_set[count, :])
    elev_bathy = elev_data_set[count, :] + bathy
    elev_bathy[elev_bathy < 0.01] = 0.0
    man = manningdg.dat.data[:]
    tau_b = np.array(1024*9.81*man*man*speed*speed / (elev_bathy)**(1./3.))
    tau_b[ elev_bathy < 0.001] = 0.0 # we have < 1mm of water
    tau_b[ tau_b < 0.0 ] = 0.0 # we had no water (shouldn't happen due to above, but just in case)
    bss.dat.data[:] = tau_b
    with CheckpointFile(output_dir + "/bss_{:05}.h5".format(count), 'w') as bss_chk:
        bss_chk.save_mesh(thetis_mesh)
        bss_chk.save_function(bss)
    bss_file.write(bss)
    count += 1

ave_speed = [] # average over speeds
max_speed = [] # max over speeds
ave_bss = [] # ave of bss calc
max_bss = [] # max of bss calc
ave_vel = [] # vector of ave u and ave v
max_vel = [] # vector of when max speed occurs
max_fs = [] # maximum wave height

for i in range(uv.dat.data.shape[0]): # loop over nodes in the Function mesh
    man = np.array(manningdg.dat.data[i])
    bathy = np.array(bathydg.dat.data[i])
    u_vel = np.array(u_data_set[:,i]) # timeseries of u, v and elevation
    v_vel = np.array(v_data_set[:,i])
    all_elev = np.array(elev_data_set[:,i])
    speed = np.sqrt(u_vel*u_vel + v_vel*v_vel)
    ave_speed.append(np.mean(speed))
    max_speed.append(np.max(speed))
    elev_bathy = elev+bathy
    power = np.sign(elev_bathy) * (np.abs(elev_bathy)) ** (1. / 3.)
    tau_b = np.array(1024*9.81*man*man*speed*speed / power)
    tau_b[ elev_bathy < 0.01] = 0.0
    tau_b[ tau_b < 0.0 ] = 0.0
    ave_bss.append(np.mean(tau_b))
    max_bss.append(np.max(tau_b))

    ave_vel.append([np.mean(u_vel), np.mean(v_vel)])
    max_vel.append([u_vel[np.argmax(speed)],v_vel[np.argmax(speed)]])

    max_fs.append(np.max(all_elev))

# We then save all the scalar temporal stats in a single hdf5 file
with CheckpointFile(output_dir + '/temporal_stats_scal.h5', "w") as chk:
    chk.save_mesh(thetis_mesh)
    avespeed = Function(P1DG, name="AveSpeed")
    avespeed.dat.data[:] = np.array(ave_speed)
    File( output_dir + '/ave_speed.pvd').write(avespeed)
    chk.save_function(avespeed)
    maxspeed = Function(P1DG, name="MaxSpeed")
    maxspeed.dat.data[:] = np.array(max_speed)
    File( output_dir + '/max_speed.pvd').write(maxspeed)
    chk.save_function(maxspeed)
    avebss = Function(P1DG, name="AveBSS")
    avebss.dat.data[:] = np.array(ave_bss)
    File( output_dir + '/average_bss.pvd').write(avebss)
    chk.save_function(avebss)
    maxbss = Function(P1DG, name="MaxBSS")
    maxbss.dat.data[:] = np.array(max_bss)
    File( output_dir + '/max_bss.pvd').write(maxbss)
    chk.save_function(maxbss, name='MaxBSS')
    maxfs = Function(P1DG, name="MaxFS")
    maxfs.dat.data[:] = np.array(max_fs)
    chk.save_function(maxfs, name='MaxFS')
    File( output_dir + '/max_fs.pvd').write(maxfs)
   
# now the vectors
with CheckpointFile(output_dir + '/temporal_stats_vec.h5', "w") as chk:
    chk.save_mesh(thetis_mesh)
    P1DG = VectorFunctionSpace(thetis_mesh, "DG", 1)
    avevel = Function(P1DG, name='ave_vel')
    avevel.dat.data[:] = np.array(ave_vel)
    File( output_dir + '/average_vel.pvd').write(avevel)
    chk.save_function(avevel, name='AveVel')
    maxvel = Function(P1DG, name='max_vel')
    maxvel.dat.data[:] = np.array(max_vel)
    File( output_dir + '/max_vel.pvd').write(maxvel)
    chk.save_function(maxvel, name='MaxVel')

