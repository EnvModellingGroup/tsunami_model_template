from thetis import *
import math
import numpy as np
import time
import sys
import slide_movement

outputdir = "output"


# Bathymetry
mesh2d = Mesh("../surface.msh")
xy = SpatialCoordinate(mesh2d)
P1_2d = FunctionSpace(mesh2d, 'CG', 1)

def depth(X):
  if (X[0] > 7000):
    depth = 400.0
  else:
    import math
    slope = math.tan(math.radians(4))
    depth = 400.0+(X[0]-7000)*slope
    if depth < 12.0:
      depth = 12.0
  return depth

bathymetry2d = Function(P1_2d, name="bathymetry")
xvector = mesh2d.coordinates.dat.data
bvector = bathymetry2d.dat.data
assert xvector.shape[0]==bvector.shape[0]
for i, (xy) in enumerate(mesh2d.coordinates.dat.data):
    bvector[i] = depth(xy)


# Set up extra simulation parameters
dt = 1                              # Time step (10 s)
t_export = 1                        # Export time (600 s)
t_end = 100                         # Run duration
h_viscosity = 1                     # Horizontal viscosity
mu_manning = 0.025                  # Manning coefficient
w_d_alpha = 10.0

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
#options.timestepper_type = 'CrankNicolson'
options.swe_timestepper_type = 'BackwardEuler'
#options.check_volume_conservation_2d = True
options.fields_to_export = ['uv_2d', 'elev_2d']
options.volume_source_2d = Function(P1_2d)
#options.use_wetting_and_drying = True
#options.wetting_and_drying_alpha = Constant(w_d_alpha)
options.horizontal_viscosity = Constant(h_viscosity)
options.use_grad_div_viscosity_term = False
options.use_grad_depth_viscosity_term = False
options.manning_drag_coefficient = Constant(mu_manning)
solver_obj.assign_initial_conditions(uv=as_vector((1e-5, 0.0)), elev = Constant(0.0))
solver_obj.bnd_functions['shallow_water'] = {
    7: {'un': 0.0},
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

