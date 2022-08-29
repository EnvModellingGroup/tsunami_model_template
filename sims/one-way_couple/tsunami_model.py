from thetis import *
import boundary_forcing
import utm
from math import sin,pi

mesh_file = '../../mesh/west_scotland_regional_hires.msh'
# things to alter (relative path to this script)
cent_lat=57.54 
cent_lon = -6.49
physID = 666
#timestepping options
dt = 6.0 # reduce this if solver does not converge
t_export = 60.0
t_start = 9500
t_end = 15.*60.*60. - t_start # 15 hours
output_directory = "output"
# model options
wetting_and_drying = True
manning_drag_coefficient = 0.025

output_dir = create_directory(output_directory)
mesh2d = Mesh(mesh_file)
print_output('Loaded mesh '+mesh2d.name)
print_output('Exporting to '+output_dir)

# bathymetry
P1 = FunctionSpace(mesh2d, "CG", 1)
# read bathymetry code
chk = DumbCheckpoint('bathymetry', mode=FILE_READ)
bathymetry2d = Function(P1)
chk.load(bathymetry2d,  name='bathymetry')
chk.close()

#read viscosity / manning boundaries code
chk = DumbCheckpoint('viscosity', mode=FILE_READ)
h_viscosity = Function(bathymetry2d.function_space(), name='viscosity')
chk.load(h_viscosity)
chk.close()

if wetting_and_drying:
    use_wetting_and_drying = True
else:
    bathymetry2d.dat.data[bathymetry2d.dat.data < 10.0] = 10.0
    use_wetting_and_drying = False

def coriolis(mesh, lat, lon):
    R = 6371e3
    Omega = 7.292e-5
    lat_r = lat * pi / 180.
    f0 = 2 * Omega * sin(lat_r)
    beta = (1 / R) * 2 * Omega * cos(lat_r)
    x = SpatialCoordinate(mesh)
    x_0, y_0, utm_zone, zone_letter = utm.from_latlon(lat, lon)
    coriolis_2d = Function(FunctionSpace(mesh, 'CG', 1), name="coriolis_2d")
    coriolis_2d.interpolate(f0 + beta * (x[1] - y_0))
    return coriolis_2d

coriolis_2d = coriolis(mesh2d, cent_lat, cent_lon)

# --- create solver ---
solverObj = solver2d.FlowSolver2d(mesh2d, bathymetry2d)
options = solverObj.options
options.no_exports = False
options.use_nonlinear_equations = True
options.simulation_export_time = t_export
options.simulation_end_time = t_end
options.output_directory = output_dir
options.check_volume_conservation_2d = True
options.fields_to_export = ['uv_2d','elev_2d']
options.fields_to_export_hdf5 = ['uv_2d', 'elev_2d']
options.element_family = "dg-dg"
options.swe_timestepper_type = 'DIRK22'
options.manning_drag_coefficient = Constant(manning_drag_coefficient)
options.coriolis_frequency = coriolis_2d
options.horizontal_viscosity = Constant(10)*h_viscosity
options.timestep = dt
options.use_wetting_and_drying = use_wetting_and_drying
if use_wetting_and_drying:
    options.use_automatic_wetting_and_drying_alpha = True
    options.wetting_and_drying_alpha_min = Constant(1.0)
    options.wetting_and_drying_alpha_max = Constant(75.0)
    #options.wetting_and_drying_alpha = Constant(wd_alpha)
options.swe_timestepper_options.solver_parameters = {
     'snes_rtol': 1e-5,
     'snes_max_it': 100,
     'snes_type': "ksponly",
     'ksp_type': 'gmres',
}

# boundary conditions
tsunami_elev = Function(FunctionSpace(mesh2d, "CG", 1), name='tsunami_elev')
solverObj.bnd_functions['shallow_water'] = {
    physID: {'elev': tsunami_elev},
    1000: {'un': 0.0},
}

# Set up usual SWE terms
solverObj.create_equations()

def update_forcings(t):
    t += t_start
    boundary_forcing.set_tsunami_field(tsunami_elev, t)

update_forcings(0.0)
solverObj.assign_initial_conditions(uv=Constant(("1e-7","0.0")), elev=Constant(0.))

# Run model
solverObj.iterate(update_forcings=update_forcings)


