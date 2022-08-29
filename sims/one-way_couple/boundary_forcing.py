import utm
import scipy.interpolate
import numpy as np

# edit me!
utm_zone=29
utm_band='U'
fluidity_times = np.loadtxt('../../data/times.txt')
boundary_coords = np.loadtxt('../../data/coordinates.txt')
boundary_data = np.loadtxt('../../data/elev.txt')


def set_tsunami_field(elev, t):
    
    # Set elev field using nearest neighbour interpolation from boundary points extracted from Fluidity tsunami simulations
    fluidity_dt = fluidity_times[1]-fluidity_times[0]
    if t%fluidity_dt == 0:
        fluidity_index = int(t/fluidity_dt)
        bd = boundary_data[fluidity_index]
        interpolator = scipy.interpolate.NearestNDInterpolator(boundary_coords, bd)
        mesh2d = elev.function_space().mesh()
        xvector = mesh2d.coordinates.dat.data
        evector = elev.dat.data
        for i,xy in enumerate(xvector):
            evector[i] = interpolator((xy[0], xy[1]))
    else:
        fluidity_index = int(t/fluidity_dt)
        boundary_data_1 = boundary_data[fluidity_index]
        boundary_data_2 = boundary_data[fluidity_index+1]
        interpolator_1 = scipy.interpolate.NearestNDInterpolator(boundary_coords, boundary_data_1)
        interpolator_2 = scipy.interpolate.NearestNDInterpolator(boundary_coords, boundary_data_2)
        mesh2d = elev.function_space().mesh()
        xvector = mesh2d.coordinates.dat.data
        evector = elev.dat.data
        for i,xy in enumerate(xvector):
            evector[i] = interpolator_1((xy[0], xy[1])) + (t-fluidity_dt*fluidity_index)*(interpolator_1((xy[0], xy[1]))-interpolator_2((xy[0], xy[1])))/fluidity_dt

