from utm import *

# path relative to the root dir of this sim. Leave as mesh/blah.msh in most cases
mesh_file = 'mesh/my_mesh.msh'
forcing_boundary = 666
utm_zone = 30
utm_band = "U"
cent_lat = 55.696
cent_lon = -1.812
start_time = 0 # linked to the forcing data. Of nothing happens in that for 3 hours, then alter this
end_time = 54000 # 15 hours
output_dir = "output"
output_time = 60 # 1 minute
