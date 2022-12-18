#!/bin/bash
# pull directory as argument
directory=$1
ncore=$2
mesh=$3

resolution=100
projection=EPSG:32630

# These are the function names in the h5 file
declare -a varname=("MaxSpeed"
                    "AveSpeed"
                    "AveBSS"
                    "MaxBSS"
                    "MaxFS"
                    )

# The English equivalent of above - *same order*, include units
declare -a names=("Max speed (m/s)"
                  "Ave speed (m/s)"
                  "Ave BSS (kgm-1s-2)"
                  "Max BSS (kgm-1s-2)"
                  "Max Wave (m)"
                 )

# loop over variables
for (( i=0; i<${#varname[@]}; i++));
do
    var=${varname[$i]}
    name=${names[$i]}
    file="${directory}/temporal_stats_scal.h5"
   	# loop over variables with counter
    echo "   Rasterising ${var}"
    # create the raster ov the vtu
    mpiexec -n ${ncore} python h5_2_raster.py --resolution ${resolution} ${file} ${mesh} temp --func ${var}
    # create a filename
    filename="${directory}/${var}".nc
    #mask it
    gdalwarp -cutline "../mesh/mask_modern.shp" -s_srs ${projection} -crop_to_cutline -of NetCDF -r bilinear  -dstnodata -9999 -overwrite temp.xyz "${filename}"
    # rename the variable to something sensible
    ncrename -v Band1,"${var}" "${filename}"		
    # change the long name variable to contain decent info
    ncatted -O -a long_name,"${var}",o,c,"${name}" "${filename}"

done

