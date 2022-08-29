from thetis import *
import math
import numpy as np
import time
import sys
import slide_movement
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import params


# Bathymetry
mesh2d = Mesh(os.path.join(os.path.pardir,os.path.pardir,params.mesh_file))
# read bathymetry code
chk = DumbCheckpoint('bathymetry', mode=FILE_READ)
bathymetry2d = Function(P1)
chk.load(bathymetry2d,  name='bathymetry')
chk.close()


#timestepping options
dt = 2 # reduce if solver does not converge
t_export = params.output_time
t_end = params.end_time
output_dir = params.output_dir
utm_zone = params.utm_zone
utm_band=params.utm_band
P1 = FunctionSpace(mesh, "CG", 1)
cent_lat = params.cent_lat
cent_lon = params.cent_lon

xy = SpatialCoordinate(mesh2d)
P1_2d = FunctionSpace(mesh2d, 'CG', 1)

#read viscosity / manning boundaries code
chk = DumbCheckpoint('viscosity', mode=FILE_READ)
h_viscosity = Function(bathymetry2d.function_space(), name='viscosity')
chk.load(h_viscosity)
chk.close()
chk = DumbCheckpoint('manning', mode=FILE_READ)
manning = Function(bathymetry2d.function_space(), name='manning')
chk.load(manning)
chk.close()


# Set up extra things to export
slide_height = Function(P1_2d, name="slide_height")         # Slide movement
slide_height_file = File(outputdir + '/slide_height.pvd')

def create_hs(t):

    hs = Function(P1_2d)
    mesh2d = hs.ufl_domain()
    xy_vector = mesh2d.coordinates.dat.data
    hs_vector = hs.dat.data
    assert xy_vector.shape[0] == hs_vector.shape[0]
    for i, (xy) in enumerate(xy_vector):
        hs_vector[i] = slide_movement.slide_height(xy,t,form="eg",profile="eg")
    return hs


# Create solver
solver_obj = solver2d.FlowSolver2d(mesh2d, bathymetry2d)
options = solver_obj.options
options.element_family = 'dg-dg'
options.polynomial_degree = 1
options.use_nonlinear_equations = True
options.output_directory = outputdir
options.timestep = dt
options.simulation_export_time = t_export
options.simulation_end_time = t_end
options.swe_timestepper_type = 'DIRK22'
options.fields_to_export = ['uv_2d', 'elev_2d']
options.fields_to_export_hdf5 = ['uv_2d', 'elev_2d']
options.volume_source_2d = Function(P1_2d)
options.use_wetting_and_drying = True
options.use_automatic_wetting_and_drying_alpha = True
options.wetting_and_drying_alpha_min = Constant(0.5)
options.wetting_and_drying_alpha_max = Constant(75.0)
options.horizontal_viscosity = h_viscosity
options.use_grad_div_viscosity_term = False
options.use_grad_depth_viscosity_term = False
options.manning_drag_coefficient = mu_manning
solver_obj.assign_initial_conditions(uv=as_vector((1e-10, 0.0)), elev = Constant(0.0))
solver_obj.bnd_functions['shallow_water'] = {
    1000: {'un': 0.0},
    #set closed boundaries to zero velocity
}
# Call creator to create the function spaces
solver_obj.create_equations()

# Update slide position and associated source term
def update_forcings(t):

    landslide_source = solver_obj.timestepper.fields.get('volume_source')
          
    hs = (create_hs(t + options.timestep) - create_hs(t))/options.timestep
    landslide_source.project(hs)

    uv, elev = solver_obj.timestepper.solution.split()
    eta = Function(P1_2d).project(elev)

    # Save extra functions to separate files every export time
    if np.mod(t,t_export) == 0:
        slide_height_file.write(slide_height.project(create_hs(t)))


solver_obj.iterate(update_forcings=update_forcings)

