from thetis import *
import numpy as np
import sys
import slide_movement
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import params

#timestepping options
dt = 2 # reduce if solver does not converge
t_export = params.output_time
t_end = params.end_time
t_start = params.start_time
output_dir = params.output_dir
utm_zone = params.utm_zone
utm_band=params.utm_band
cent_lat = params.cent_lat
cent_lon = params.cent_lon

# read bathymetry code
chk = CheckpointFile('bathymetry', 'r')
mesh2d = chk.load_mesh()
bathymetry2d = chk.load_function(mesh2d,'bathymetry')
chk.close()

#read viscosity / manning boundaries code
chk = CheckpointFile('viscosity', 'r')
mesh2d = chk.load_mesh()
h_viscosity = chk.load_function(mesh2d,'viscosity')
chk.close()
chk = CheckpointFile('manning', 'r')
mesh2d = chk.load_mesh()
manning = chk.load_function(mesh2d, 'manning')
chk.close()

xy = SpatialCoordinate(mesh2d)
P1_2d = FunctionSpace(mesh2d, 'CG', 1)

# function to set up the Coriolis force
# Depends on a "central" lat/lon point in
# your mesh
def coriolis(mesh, lat, lon):
    R = 6371e3
    Omega = 7.292e-5
    lat_r = lat * pi / 180.
    f0 = 2 * Omega * sin(lat_r)
    beta = (1 / R) * 2 * Omega * cos(lat_r)
    x = SpatialCoordinate(mesh)
    x_0, y_0, utm_zone, zone_letter = params.from_latlon(lat, lon)
    coriolis_2d = Function(FunctionSpace(mesh, 'CG', 1), name="coriolis_2d")
    coriolis_2d.interpolate(f0 + beta * (x[1] - y_0))

    return coriolis_2d

# Set up extra things to export
slide_height = Function(P1_2d, name="slide_height")         # Slide movement
slide_height_file = File(outputdir + '/slide_height.pvd')

#account for Coriolis code - mesh, centre lat, centre lon
coriolis_2d = coriolis(mesh, cent_lat, cent_lon)

def create_hs(t):
    hs = Function(P1_2d)
    mesh2d = hs.ufl_domain()
    xy_vector = mesh2d.coordinates.dat.data
    hs_vector = hs.dat.data
    assert xy_vector.shape[0] == hs_vector.shape[0]
    for i, xy in enumerate(xy_vector):
        hs_vector[i] = set_slide_height(xy,t,form=params.FORM,profile=params.PROFILE)
    return hs


# --- create solver ---
solverObj = solver2d.FlowSolver2d(mesh, bathymetry2d)
options = solverObj.options
options.use_nonlinear_equations = True
options.simulation_export_time = t_export
options.simulation_end_time = t_end
options.output_directory = output_dir
options.check_volume_conservation_2d = True
options.fields_to_export = ['uv_2d', 'elev_2d']
options.fields_to_export_hdf5 = ['uv_2d', 'elev_2d']
options.manning_drag_coefficient = manning #the manning function we created in initialisation & loaded above
options.horizontal_viscosity = h_viscosity #the viscosity 'cushion' we created in initialisation & loaded above
options.coriolis_frequency = coriolis_2d
options.timestep = dt
options.volume_source_2d = Function(P1_2d)
options.use_automatic_wetting_and_drying_alpha = True
options.wetting_and_drying_alpha_min = Constant(0.5)
options.wetting_and_drying_alpha_max = Constant(75.0)
options.use_wetting_and_drying = True
options.element_family = "dg-dg"
options.swe_timestepper_type = 'DIRK22'


solverObj.bnd_functions['shallow_water'] = {
    1000: {'un': 0.0},
    #set closed boundaries to zero velocity
}

solverObj.assign_initial_conditions(uv=as_vector((1e-10, 0.0)), elev = Constant(0.0))

# Call creator to create the function spaces
solverObj.create_equations()

# Update slide position and associated source term
def update_forcings(t):

    landslide_source = solverObj.timestepper.fields.get('volume_source')
          
    hs = (create_hs(t + options.timestep) - create_hs(t))/options.timestep
    landslide_source.project(hs)

    uv, elev = solverObj.timestepper.solution.split()
    eta = Function(P1_2d).project(elev)

    # Save extra functions to separate files every export time
    if np.mod(t,t_export) == 0:
        slide_height_file.write(slide_height.project(create_hs(t)))


solverObj.iterate(update_forcings=update_forcings)
