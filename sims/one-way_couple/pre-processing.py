import hrds
from firedrake import *
from firedrake.petsc import PETSc


PETSc.Sys.Print('setting up mesh across %d processes' % COMM_WORLD.size)
manning_drag = 0.025
mesh_file = '../../mesh/west_scotland_regional_hires.msh'
bathymetry_files = ["../../data/palaeo_bathytopo_50m_utm29.tif","../../data/palaeo_lidar_data_nSkye_final_utm.tif", "../../data/palaeo_lidar_data_nanCorr_final_utm.tif", "../../data/palaeo_lidar_data_gorten_final_utm.tif"]
forcing_boundary = 666
distances = [100.0, 100.0, 100.0]
viscosity = 1.0
factor = 1000.
length_scale = 6.0e3
mesh2d = Mesh(mesh_file)
PETSc.Sys.Print("mesh loaded")

def get_bathymetry(bathymetry_files, mesh2d, distances, minimum_depths=None, negative_depth=True):
    P1_2d = FunctionSpace(mesh2d, 'CG', 1)
    bathymetry_2d = Function(P1_2d, name="bathymetry")
    xvector = mesh2d.coordinates.dat.data
    bvector = bathymetry_2d.dat.data
    assert xvector.shape[0]==bvector.shape[0]

    bathy = hrds.HRDS(bathymetry_files[0], 
             rasters=bathymetry_files[1:], 
             distances=distances,saveBuffers=False)
    bathy.set_bands()
    for i, (xy) in enumerate(mesh2d.coordinates.dat.data):
        bvector[i] = -1.0 * bathy.get_val(xy)
    
    return bathymetry_2d

def smoothen_bathymetry(bathymetry2d): # smoothing bathymetry
    v = TestFunction(bathymetry2d.function_space())
    massb = assemble(v * bathymetry2d *dx)
    massl = assemble(v*dx)
    with massl.dat.vec as ml, massb.dat.vec as mb, bathymetry2d.dat.vec as sb:
        ml.reciprocal()
        sb.pointwiseMult(ml, mb)

# first deal with bathymetry
chk = DumbCheckpoint('bathymetry', mode=FILE_CREATE)
bathymetry2d = get_bathymetry(bathymetry_files, mesh2d, distances=distances)
smoothen_bathymetry(bathymetry2d)
chk.store(bathymetry2d, name='bathymetry')

# now create distance from boundary function
# typical length scale
L = 1e3
V = FunctionSpace(mesh2d, 'CG', 1)
# Calculate distance to open boundary
PETSc.Sys.Print('Calculate distance for viscosity')
bcs = [DirichletBC(V, 0.0, forcing_boundary)] #make sure this matches physicalID of open boundaries
v = TestFunction(V)
u = Function(V)
solver_parameters = {
    'ksp_rtol': 1e-3,
}
# Before we solve the Eikonal equation, let's solve a Laplace equation to
# generate an initial guess
F = L**2*(inner(grad(u), grad(v))) * dx - v * dx
PETSc.Sys.Print('solve 1')
solve(F == 0, u, bcs, solver_parameters=solver_parameters)


solver_parameters = {
    'snes_type': 'newtonls',
    'snes_monitor': None,
    'ksp_rtol': 1e-4,
    'ksp_type': 'preonly',
    "ksp_max_it": 2000,
    "ksp_converged_reason": None,
    'pc_type': 'lu',
    'pc_factor_mat_solver_packages': 'mumps',
    "ksp_view": None,
    "ksp_monitor_true_residual": None,
}


# epss values set the accuracy (in meters) of the final "distance to boundary" function. To make
# more accurate add in extra iterations, eg, 500., 250., etc. This may result in the solver not
# converging.
epss = [100000., 10000., 7500., 2500., 1500., 1000.]
# solve Eikonal equations
for i, eps in enumerate(epss):
    PETSc.Sys.Print('Solving Eikonal with eps == ' + str(float(eps)))
    F = inner(sqrt(inner(grad(u), grad(u))), v) * dx - v * dx + eps*inner(grad(u), grad(v)) * dx
    solve(F == 0, u, bcs, solver_parameters=solver_parameters)

chk = DumbCheckpoint('dist', mode=FILE_CREATE)
dist = Function(V, name='dist')
dist.interpolate(u)
chk.store(dist, name='dist')
File('dist.pvd').write(u)

# create a viscosity buffer
#create boundary of increased viscosity
chk = DumbCheckpoint('viscosity', mode=FILE_CREATE)
h_viscosity = Function(V, name='viscosity')
h_viscosity.interpolate(Max(viscosity, (viscosity*factor) * (viscosity - u / length_scale)))
chk.store(h_viscosity, name='viscosity')
File('viscosity.pvd').write(h_viscosity)

