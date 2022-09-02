import numpy as np
from thetis import *
import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import params

# which mesh to load?
mesh2d = Mesh(os.path.join(os.path.pardir,os.path.pardir,params.mesh_file))

# where should the output of this analysis go
output_dir = 'analysis'
create_directory(output_dir)

# where did thetis put the hdf5 files?
file_location = params.output_dir + '/hdf5/' #location of the Elevation2d output files

# How long does your simulations run for (s)
t_end = params.end_time #40 days (i.e. 30 days of analysis)
# how often are exports produced in the main run?
t_export = params.output_time
# which is the start file?
start_file = int(params.spin_up / t_export)

# You shouldn't need to edit below here
#========================================
xvector = mesh.coordinates.dat.data

t_n = int(t_end/t_export + 1)
thetis_times = t_export*np.arange(t_n) + t_export
P1 = FunctionSpace(mesh, "CG", 1)
P2 = FunctionSpace(mesh, "DG", 1)
elev = Function(P2, name='elev_2d')
elev_data_set = np.empty((t_n, elev.dat.data.shape[0]))
for i in range(start_file,t_n):
    print('Reading h5 files. Time ',i,i*t_export)
    checkpoint_file = DumbCheckpoint(file_location + '/Elevation2d_{:05}'.format(i),    mode=FILE_READ)
    checkpoint_file.load(elev)
    checkpoint_file.close()
    elev_data_set[i, :] = elev.dat.data[:]

detector_maxfs = []
detector_minfs = []
for i in range(elev.dat.data.shape[0]):
    thetis_elev = elev_data_set[:,i]
    
    detector_maxfs.append(max(thetis_elev[int(start_file):]))
    detector_minfs.append(min(thetis_elev[int(start_file):]))

# sort out the min, max and tidal range - save as h5 to rasterise.
minfs = Function(P2, name="MinFS")
minfs.dat.data[:] = np.array(detector_minfs)
chk = DumbCheckpoint(output_dir + '/min_fs', mode=FILE_CREATE)
chk.store(minfs, name='MinFS')
File( output_dir + '/min_fs.pvd').write(minfs)
maxfs = Function(P2, name="MaxFS")
maxfs.dat.data[:] = np.array(detector_maxfs)
chk = DumbCheckpoint(output_dir + '/max_fs', mode=FILE_CREATE)
chk.store(maxfs, name='MaxFS')
File( output_dir + '/max_fs.pvd').write(maxfs)

