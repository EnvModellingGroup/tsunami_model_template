#!/bin/bash
# pull directory as argument
directory=$1
ncore=$2
mesh=$3

resolution=5000
projection=EPSG:32630

# filenames to rasterise, without .h5 extension
declare -a varname=("max_speed"
                    "max_bss"
                    "max_fs"
                    )

# The English equivalent of above - *same order*, include units
declare -a names=("Max speed (m/s)"
                  "Max BSS (kgm-1s-2)"
                  "Max Wave (m)"
                 )

# loop over variables
for (( i=0; i<${#varname[@]}; i++));
do
    var=${varname[$i]}
    name=${names[$i]}
    echo ${var}
    file="${directory}/${var}"
    echo "Dealing with ${file}"
   	# loop over variables with counter
    echo "   Rasterising ${var}"
    # create the raster ov the vtu
    mpiexec -n ${ncore} python h5_2_raster.py --resolution ${resolution} ${file} ${mesh} temp
    # create a filename
    filename="${directory}/${var}".nc
    #mask it
    gdalwarp -cutline "../mesh/mask.shp" -s_srs ${projection} -crop_to_cutline -of NetCDF -r bilinear  -dstnodata -9999 -overwrite temp.xyz "${filename}"
    # rename the variable to something sensible
    ncrename -v Band1,"${var}" "${filename}"		
    # change the long name variable to contain decent info
    ncatted -O -a long_name,"${var}",o,c,"${name}" "${filename}"

done

