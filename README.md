
This is my work from summer 2016. The goal is to add proper XDMF support to VisIt for use in the Chimera collaboration. It builds on the work Jake Rosa did during summer 2014. In this installation two paired large files managed by git lfs that can also be found on the UTK newton server in /lustre/projects/astro/chimera/lentz/. The files are chimera_00774_grid_1_01.h5, and d96-2d-sn160-00774.silo. The SILO file was created with a script from the HDF5 file and displays properly in VisIt. We seek to eliminate the need for this script by using XDMF. With this goal in mind, two python scripts now exist within this repository: write_xml.py and reducer.py. 

write_xml.py
============
Write_xml.py takes input HDF files from chimera and makes XDMF to read them. It also optionally reduces them down to a size Bellarophon can handle by copying over only the hydro and abundance scalars. Radiation quantities are computed optionally and also added. If a reduction step is not taken but the option to use computed quantities is selected, the computed quantities will be placed in another HDF file with the same name as the input file only without 

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
Equally valid to VisIt is the __.xmdf__ extension. I prefer shorter extensions usually but this is more descriptive so one can use can use this extension instead of __.xmf__ by adding the `--xdmf` flag. Additionally, if one would rather no have the words _step_ and _grid_ in the file name, these can be eliminated by using the `--short` flag. Calling the script with both these options will create a file with the same contents as before but now will the file name `chimera-1-003800000.xdmf`

**Options:**

* __`--quiet`__   
This option was created so the script can run "headless" and will only output to the stderr data stream in the event of a problem. If a human calls this script, it's probably better not to use this option as if something breaks or hangs slightly or unpredictably, one won't have much idea where or why.

* __`--slices int`__   
This option will limit the number of wedges used for diagnostic reasons. For example on the 003800000 time step there are 30 slices and each is ~1.8GB so to have them all in a directory requires 54GB of space. This might be a bit much to download so if you only download the first 5 or so you can limit the script to use just those by calling `$python write_xml.py chimera_003800000_grid_1_01.h5 --slices 5`. The output XDMF file when opened with VisIt will just display those wedges instead of the entire sphere. 

* __`--norepeat`__   
Disregard this option. It's internally needed to prevent the script from infinitely calling itself when reloading modules on Rhea if something breaks really badly. If the wrong modules are loaded on Rhea or Titan the script will automatically spawn the a sub shell, load the right modules, and then call itself. The problem is if it doesn't work and this flag isn't used it will just continue to do this forever so this option makes it so that it won't spawn any sub-shells and ends the madness. 

* __`--threads int`__   
For the auxiliary/computed variables the script runs with up to 16 cores. Unfortunately I currently have no way of distinguishing real from logical cores and if the script runs with more cores than there are real cores it will actually run dramatically slower. This option simply provides some more manual control over how the script runs. On Rhea there were some permission changes recently that broke the multiprocessing module so the this works as a manual override by just setting the threads=1 with `--threads 1`.

* __`--prefix str`__ or __`-p str`__   
Suppose you have the input filename as `foo_0012345_grid_1_01.h5` and use the command `$python write_xml.py foo_0012345_grid_1_01.h5 -s`. Ordinarily one would get the file `foo-0012345-1.xmf` but if one instead called: `$python write_xml.py foo_0012345_grid_1_01.h5 -s -prefix bar` one would instead get `bar-0012345-1.xmf`.

* __`--disable str [str]`__   
Will skip writing XDMF for any of the following specified sections:
    * `abundance`
    * `hydro`
    * `radiation` 
They can also be combined like so:
    `$python write_xml.py foo_12345_grid_1_01.h5 --disable abundance radiation`
This would only write xdmf for the hydro variables.