
This is my work from summer 2016. The goal is to add proper XDMF support to VisIt for use in the Chimera collaboration. It builds on the work Jake Rosa did during summer 2014. In this installation two paired large files managed by git lfs that can also be found on the UTK newton server in /lustre/projects/astro/chimera/lentz/. The files are chimera_00774_grid_1_01.h5, and d96-2d-sn160-00774.silo. The SILO file was created with a script from the HDF5 file and displays properly in VisIt. We seek to eliminate the need for this script by using XDMF. With this goal in mind, two python scripts now exist within this repository: write_xml.py and reducer.py. 

write_xml.py
============
Write_xml.py takes input HDF files from chimera and makes XDMF to read them. It also optionally reduces them down to a size Bellarophon can handle by copying over only the hydro and abundance scalars. Radiation quanties are computed optionally and also added. If a reduction step is not taken but the option to use computed quanties is selected, the computed quanties will be placed in another HDF file with the same name as the input file only without 

###Usage
    usage: write_xml.py [-h] [--threads int] [--slices [int]] [--prefix [str]]
                        [--repeat] [--quiet] [--short] [--norepeat]
                        [--disable str [str ...]] [--xdmf] [--directory [str]]
                        [--auxiliary]
                        foo.h5 [foo.h5 ...]
    
    Generate XDMF files from Chimera hdf5 files
    
    positional arguments:
      foo.h5                hdf5 files to process (1 or more args)
    
    optional arguments:
      -h, --help            show this help message and exit
      --threads int         specify number of threads
      --slices [int]        number of slices to use
      --prefix [str], -p [str]
                            specify the xmf file prefix
      --repeat, -r          use the first wedge for all slices
      --quiet, -q           only display error messages (default full debug
                            messages)
      --short, -s           use shorter filenaming convention
      --norepeat            debug variable for infinite recursive execution
                            escaping
      --disable str [str ...]
                            disable output of abundances, hydro, or radiation
                            components
      --xdmf                use .xdmf extension instead of default .xmf
      --directory [str]     Output xdmf in dirctory specified instead of next 
                            to hdf files
      --auxiliary, -a       Write auxiliary computed (derivative) values like
                            luminosity to a companion file
An example call might be like:

    $python write_xml.py chimera_003800000_grid_1_01.h5

This would just create a file titled `chimera_grid-1_step-003800000.xmf`
Equally valid to VisIt is the __.xmdf__ extension. I prefer shorter extensions usually but this is more descriptive so one can use can use this extension instead of __.xmf__ by adding the `--xdmf` flag. Additionaly if one would rather no have the words _step_ and _grid_ in the file name, these can be eliminated by using the `--short` flag. Calling the script with both these options will create a file with the same contents as before but now will the file name `chimera-1-003800000.xdmf`

Other options:

**`--quiet`**
   This option was created 