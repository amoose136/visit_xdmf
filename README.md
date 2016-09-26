
This is my work from summer 2016. The goal is to add proper XDMF support to VisIt for use in the Chimera collaboration. It builds on the work Jake Rosa did during summer 2014. In this installation is a pair of large files managed by Git LFS that can also be found on the UTK newton server in /lustre/projects/astro/chimera/lentz/. The files are chimera_00774_grid_1_01.h5, and d96-2d-sn160-00774.silo. The SILO file was created with a script from the HDF5 file and displays properly in VisIt. We seek to eliminate the need for this script by using XDMF. With this goal in mind, two python scripts now exist within this repository: write_xml.py and reducer.py. 

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

####Options
* __`--short`__ or __`-s`__   
If one would rather not have the words _step_ and _grid_ in the file name, these can be eliminated by using the `--short` flag.
Calling `$python write_xml.py foo_0123_grid_1_01.h5 -s` will yield a xdmf file name titled `foo-1-0123.xmf` instead of `foo_grid-1_step-0123.xmf` as would be default.

* __`--quiet`__ or __`-q`__  
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
    + `abundance`
    + `hydro`
    + `radiation`

    They can also be combined like so:
    `$python write_xml.py foo_12345_grid_1_01.h5 --disable abundance radiation`
    This would only write xdmf for the hydro variables.

* __`--repeat`__ or __`-r`__   
If this option is invoked the data from the first wedge will be used for all wedges. This is useful for diagnostic reasons.

* __`--auxiliary`__ or __`-a`__   
If this option is invoked the script will also create an auxiliary HDF file containing the computed quantities.

* __`--directory [str]`__   
By default the script will output the XDMF file in the same directory as the input HDF file. This allows you to override by providing a path relative to the terminal location at the time the script is called. Example:
    `$python write_xml.py ../foo_12345_grid_1_01.h5 --directory .` 
will output the XDMF file to the current directory instead of `..`.

* __`--xdmf`__   
Equally valid to VisIt is the __.xmdf__ extension. I prefer shorter extensions usually but this is more descriptive so one can use this extension instead of __.xmf__ by adding the `--xdmf` flag.


####Examples
With shorter style filename option and also the xdmf flag (headless):
**shell:** `$python write_xml.py foo_123_grid_1_01.h5 -qs --xdmf`
**output:** `foo-1-123.xdmf` is created
___
With computed auxiliary scalars like E_RMS_[et,et-bar,mt,mt-bar] and Luminosity_[et,et-bar,mt,mt-bar] while also using the shorter naming convention (headless): 
**shell:** `$python write_xml.py foo_007_grid_2_01.h5 -qas`
**output:** `foo-2-007.xmf`,`foo_007_grid_2_aux.h5` are created
___
Using only the first slice (headless):
**shell:** `$python write_xml.py foo_123_grid_1_01.h5 -q --slices 1`
**output:** `foo_grid-1_step-123.xmf` is created
___
Using the data from first wedge for the entire grid (headless):
**shell:** `$python write_xml.py foo_123_grid_1_01.h5 -qr`
**output:** `foo_grid-1_step-123.xmf` is created
___
Using the data from first wedge for the entire grid:
**shell:** `$python write_xml.py foo_123_grid_1_01.h5 -r`
**output:**
STDOUT: 
```
Running with lxml.etree
Running with single thread
--- foo_grid-1_step-123.xmf created in 0.657558917999 seconds ---
```

`foo_grid-1_step-123.xmf` is created
___
Using the data from the first wedge for the entire grid and compute auxiliary scalars like E_RMS_[et,et-bar,mt,mt-bar] and Luminosity_[et,et-bar,mt,mt-bar] and use shorter naming convention:
**shell:** `$python write_xml.py foo_123_grid_1_01.h5 -rsa`
**output:** 
STDOUT: 
```
Running with lxml.etree
Running with single thread
Creating derived values
Computing E_RMS_[1..4] for slice 1 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 2 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 3 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 4 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 5 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 6 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 7 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 8 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 9 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 10 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 11 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 12 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 13 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 14 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 15 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 16 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 17 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 18 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 19 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 20 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 21 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 22 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 23 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 24 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 25 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 26 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 27 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 28 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 29 from foo_123_grid_1_01.h5
Computing E_RMS_[1..4] for slice 30 from foo_123_grid_1_01.h5
Writing E_RMS_0 out to file
Writing E_RMS_1 out to file
Writing E_RMS_2 out to file
Writing E_RMS_3 out to file
Computing luminosities
Computing luminosity for species 0:
    On slice 1 of 30 from foo_123_grid_1_01.h5
    On slice 2 of 30 from foo_123_grid_1_01.h5
    On slice 3 of 30 from foo_123_grid_1_01.h5
    On slice 4 of 30 from foo_123_grid_1_01.h5
    On slice 5 of 30 from foo_123_grid_1_01.h5
    On slice 6 of 30 from foo_123_grid_1_01.h5
    On slice 7 of 30 from foo_123_grid_1_01.h5
    On slice 8 of 30 from foo_123_grid_1_01.h5
    On slice 9 of 30 from foo_123_grid_1_01.h5
    On slice 10 of 30 from foo_123_grid_1_01.h5
    On slice 11 of 30 from foo_123_grid_1_01.h5
    On slice 12 of 30 from foo_123_grid_1_01.h5
    On slice 13 of 30 from foo_123_grid_1_01.h5
    On slice 14 of 30 from foo_123_grid_1_01.h5
    On slice 15 of 30 from foo_123_grid_1_01.h5
    On slice 16 of 30 from foo_123_grid_1_01.h5
    On slice 17 of 30 from foo_123_grid_1_01.h5
    On slice 18 of 30 from foo_123_grid_1_01.h5
    On slice 19 of 30 from foo_123_grid_1_01.h5
    On slice 20 of 30 from foo_123_grid_1_01.h5
    On slice 21 of 30 from foo_123_grid_1_01.h5
    On slice 22 of 30 from foo_123_grid_1_01.h5
    On slice 23 of 30 from foo_123_grid_1_01.h5
    On slice 24 of 30 from foo_123_grid_1_01.h5
    On slice 25 of 30 from foo_123_grid_1_01.h5
    On slice 26 of 30 from foo_123_grid_1_01.h5
    On slice 27 of 30 from foo_123_grid_1_01.h5
    On slice 28 of 30 from foo_123_grid_1_01.h5
    On slice 29 of 30 from foo_123_grid_1_01.h5
    On slice 30 of 30 from foo_123_grid_1_01.h5
Computing luminosity for species 1:
    On slice 1 of 30 from foo_123_grid_1_01.h5
    On slice 2 of 30 from foo_123_grid_1_01.h5
    On slice 3 of 30 from foo_123_grid_1_01.h5
    On slice 4 of 30 from foo_123_grid_1_01.h5
    On slice 5 of 30 from foo_123_grid_1_01.h5
    On slice 6 of 30 from foo_123_grid_1_01.h5
    On slice 7 of 30 from foo_123_grid_1_01.h5
    On slice 8 of 30 from foo_123_grid_1_01.h5
    On slice 9 of 30 from foo_123_grid_1_01.h5
    On slice 10 of 30 from foo_123_grid_1_01.h5
    On slice 11 of 30 from foo_123_grid_1_01.h5
    On slice 12 of 30 from foo_123_grid_1_01.h5
    On slice 13 of 30 from foo_123_grid_1_01.h5
    On slice 14 of 30 from foo_123_grid_1_01.h5
    On slice 15 of 30 from foo_123_grid_1_01.h5
    On slice 16 of 30 from foo_123_grid_1_01.h5
    On slice 17 of 30 from foo_123_grid_1_01.h5
    On slice 18 of 30 from foo_123_grid_1_01.h5
    On slice 19 of 30 from foo_123_grid_1_01.h5
    On slice 20 of 30 from foo_123_grid_1_01.h5
    On slice 21 of 30 from foo_123_grid_1_01.h5
    On slice 22 of 30 from foo_123_grid_1_01.h5
    On slice 23 of 30 from foo_123_grid_1_01.h5
    On slice 24 of 30 from foo_123_grid_1_01.h5
    On slice 25 of 30 from foo_123_grid_1_01.h5
    On slice 26 of 30 from foo_123_grid_1_01.h5
    On slice 27 of 30 from foo_123_grid_1_01.h5
    On slice 28 of 30 from foo_123_grid_1_01.h5
    On slice 29 of 30 from foo_123_grid_1_01.h5
    On slice 30 of 30 from foo_123_grid_1_01.h5
Computing luminosity for species 2:
    On slice 1 of 30 from foo_123_grid_1_01.h5
    On slice 2 of 30 from foo_123_grid_1_01.h5
    On slice 3 of 30 from foo_123_grid_1_01.h5
    On slice 4 of 30 from foo_123_grid_1_01.h5
    On slice 5 of 30 from foo_123_grid_1_01.h5
    On slice 6 of 30 from foo_123_grid_1_01.h5
    On slice 7 of 30 from foo_123_grid_1_01.h5
    On slice 8 of 30 from foo_123_grid_1_01.h5
    On slice 9 of 30 from foo_123_grid_1_01.h5
    On slice 10 of 30 from foo_123_grid_1_01.h5
    On slice 11 of 30 from foo_123_grid_1_01.h5
    On slice 12 of 30 from foo_123_grid_1_01.h5
    On slice 13 of 30 from foo_123_grid_1_01.h5
    On slice 14 of 30 from foo_123_grid_1_01.h5
    On slice 15 of 30 from foo_123_grid_1_01.h5
    On slice 16 of 30 from foo_123_grid_1_01.h5
    On slice 17 of 30 from foo_123_grid_1_01.h5
    On slice 18 of 30 from foo_123_grid_1_01.h5
    On slice 19 of 30 from foo_123_grid_1_01.h5
    On slice 20 of 30 from foo_123_grid_1_01.h5
    On slice 21 of 30 from foo_123_grid_1_01.h5
    On slice 22 of 30 from foo_123_grid_1_01.h5
    On slice 23 of 30 from foo_123_grid_1_01.h5
    On slice 24 of 30 from foo_123_grid_1_01.h5
    On slice 25 of 30 from foo_123_grid_1_01.h5
    On slice 26 of 30 from foo_123_grid_1_01.h5
    On slice 27 of 30 from foo_123_grid_1_01.h5
    On slice 28 of 30 from foo_123_grid_1_01.h5
    On slice 29 of 30 from foo_123_grid_1_01.h5
    On slice 30 of 30 from foo_123_grid_1_01.h5
Computing luminosity for species 3:
    On slice 1 of 30 from foo_123_grid_1_01.h5
    On slice 2 of 30 from foo_123_grid_1_01.h5
    On slice 3 of 30 from foo_123_grid_1_01.h5
    On slice 4 of 30 from foo_123_grid_1_01.h5
    On slice 5 of 30 from foo_123_grid_1_01.h5
    On slice 6 of 30 from foo_123_grid_1_01.h5
    On slice 7 of 30 from foo_123_grid_1_01.h5
    On slice 8 of 30 from foo_123_grid_1_01.h5
    On slice 9 of 30 from foo_123_grid_1_01.h5
    On slice 10 of 30 from foo_123_grid_1_01.h5
    On slice 11 of 30 from foo_123_grid_1_01.h5
    On slice 12 of 30 from foo_123_grid_1_01.h5
    On slice 13 of 30 from foo_123_grid_1_01.h5
    On slice 14 of 30 from foo_123_grid_1_01.h5
    On slice 15 of 30 from foo_123_grid_1_01.h5
    On slice 16 of 30 from foo_123_grid_1_01.h5
    On slice 17 of 30 from foo_123_grid_1_01.h5
    On slice 18 of 30 from foo_123_grid_1_01.h5
    On slice 19 of 30 from foo_123_grid_1_01.h5
    On slice 20 of 30 from foo_123_grid_1_01.h5
    On slice 21 of 30 from foo_123_grid_1_01.h5
    On slice 22 of 30 from foo_123_grid_1_01.h5
    On slice 23 of 30 from foo_123_grid_1_01.h5
    On slice 24 of 30 from foo_123_grid_1_01.h5
    On slice 25 of 30 from foo_123_grid_1_01.h5
    On slice 26 of 30 from foo_123_grid_1_01.h5
    On slice 27 of 30 from foo_123_grid_1_01.h5
    On slice 28 of 30 from foo_123_grid_1_01.h5
    On slice 29 of 30 from foo_123_grid_1_01.h5
    On slice 30 of 30 from foo_123_grid_1_01.h5
########################################
Writing luminosity species 0 to auxilary hdf file
Writing luminosity species 1 to auxilary hdf file
Writing luminosity species 2 to auxilary hdf file
Writing luminosity species 3 to auxilary hdf file
########################################
Creating On_grid_mask
Writing On_grid_mask to auxillary file
--- foo-1-123.xmf created in 161.676426172 seconds ---
```
`foo_123_grid_1_aux.h5` is created.
`foo-1-123.xmf` is created. (note:`foo-1-123.xmf` knows about both `foo_123_grid_1_aux.h5` and `foo_123_grid_1_`[00..N_wedges]`.h5`)
 
___
Once these files are created one needs only the the xdmf file in VisIt and it should behave almost as drop-in replacement for Silo. 

# Reducer.py
TBA 
(The short story is that it's almost the same in invocation as write_xml.py but it works only on 3d files and slices through the x,y, and z axes to create 3 2D HDF files and an xdmf file to link them together nicely. It also makes the auxiliary scalars and if you're just going to slice the data about one of these plane in VisIt anyway, this is much much faster.)
