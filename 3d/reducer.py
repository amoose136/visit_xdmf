#!/usr/bin/python
#  ______   ______  _____    _    _   ______  ______  ______  
# | |  | \ | |     | | \ \  | |  | | | |     | |     | |  | \ 
# | |__| | | |---- | |  | | | |  | | | |     | |---- | |__| | 
# |_|  \_\ |_|____ |_|_/_/  \_|__|_| |_|____ |_|____ |_|  \_\ 
# 		A script to reduce 3d datasets to 2d               
from __future__ import print_function
# For diagnostics, time the execution of the code
import time
start_time = time.time()
# Import all the things robustly
############################################################################################################################################################################################
import sys, os, math, socket, argparse, re
import subprocess as sp
import numpy as np
from pdb import set_trace as br #For debugging I prefer the c style "break" nomenclature to "trace"
import multiprocessing as mp #For parallel speedup in derivative values 
#define an error printing function for error reporting to terminal STD error IO stream
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
#define a standard printing function that only functions if there is no silence flag on script invocation
def qprint(*arg,**kwargs):
	if not args.quiet:
		print(*arg,**kwargs)
# For ORNL
if socket.gethostname()[:4]=='rhea':
	sys.path.append('/lustre/atlas/proj-shared/ast109/amos/lib/python2.7/site-packages')
	sys.path.append('/sw/redhat6/visit/current/linux-x86_64/lib/site-packages/')
if socket.gethostname()[:5]=='titan':
	sys.path.append('/lustre/atlas/proj-shared/ast109/amos/lib/python2.7/site-packages')
# For my laptop
elif socket.gethostname()=='Lycoris':
	sys.path.append('/Applications/VisIt.app/Contents/Resources/2.10.2/darwin-x86_64/lib/site-packages')

from joblib import Parallel, delayed
import six
if __name__ == '__main__':
	# Contstruct Parser:
	parser = argparse.ArgumentParser(description="Reduce 3d chimera files to 2d files")
	# parser.files is a list of 1 or more h5 files that will have xdmf files generated for them
	parser.add_argument('files',metavar='foo.h5',type=str,nargs='+',help='hdf5 files to process (1 or more args)')
	parser.add_argument('--extent','-e',dest='dimensions',metavar='int',action='store',type=int,nargs=1, help='Set max R cell to this number (by grid cell number not spatial dimensions)')
	parser.add_argument('--prefix','-p',dest='prefix',metavar='str',action='store',type=str,nargs='?', help='specify the xmf file prefix')
	parser.add_argument('--repeat','-r',dest='repeat',action='store_const',const=True, help='use the first wedge for all slices')
	parser.add_argument('--quiet','-q',dest='quiet',action='store_const',const=True, help='only display error messages (default full debug messages)')
	parser.add_argument('--short','-s',dest='shortfilename',action='store_const',const=True, help='use shorter filenaming convention')
	parser.add_argument('--norepeat',dest='norepeat',action='store_const',const=True, help='debug variable for infinite recursive execution escaping')
	parser.add_argument('--disable',dest='disable',action='store',metavar='str',type=str, nargs='+', help='disable output of abundances, hydro, or radiation components')
	parser.add_argument('--xdmf',dest='xdmf',action='store_const',const=True, help='use .xdmf extension instead of default .xmf')
	parser.add_argument('--directory',dest='dir',metavar='str', const='.',action='store',type=str,nargs='?',help='Output xdmf in dirctory specified instead of next to hdf files')
	parser.add_argument('--auxiliary','-a',dest='aux',action='store_const',const=True, help='Write auxiliary computed (derivative) values like luminosity to a companion file')
	args=parser.parse_args()
	#End Parser construction
	#####################################################################################################################################################################################################
	#This next bit is specific to ORNL. If h5py import fails it switches environments and reloads this script
	try:
		#Most likly to fail part.
		import h5py
	# Try and correct the h5import error by launching subprocess that calls this script again after loading proper modules on rhea
	except ImportError:
		try:
			qprint("Trying to run under reloaded modules")
			try:
				if not args.norepeat:
					if socket.gethostname()[:4]=='rhea':
						sp.call(['bash -cl "cd '+os.getcwd()+'; module unload PE-intel python;module load PE-gnu python python_h5py;python '+(' '.join(sys.argv))+' --norepeat"'],shell=True)
					if socket.gethostname()[:5]=='titan':
						sp.call(['bash -cl "cd '+os.getcwd()+'; module load python python_h5py;python '+(' '.join(sys.argv))+' --norepeat"'],shell=True)
				else:
					raise ValueError('aw crap, already reloaded and still failed to import something')
			except:
				#redo the offending call so the error can display
				sp.call(["module unload PE-intel python;module load PE-gnu python python_h5py"],shell=True)
				eprint("Could not import modules")
				raise ValueError('aw crap')
			qprint("Finished")
		except:
			eprint("Fatal error: could not import h5py or reload modules to make it possible. h5 reading and writing is impossible without h5py.")
		sys.exit()
	#Robustly import an xml writer/parser
	try:
		from lxml import etree as et
		qprint("Running with lxml.etree")
	except ImportError:
		try:
			# Python 2.5
			import xml.etree.cElementTree as et
			import xml.dom.minidom as md
			qprint("Running with cElementTree on Python 2.5+")
		except ImportError:
			try:
			# Python 2.5
				import xml.etree.ElementTree as et
				qprint("Running with ElementTree on Python 2.5+")
			except ImportError:
				try:
					# normal cElementTree install
					import cElementTree as et
					qprint("running with cElementTree")
				except ImportError:
					try:
						# normal ElementTree install
						import elementtree.ElementTree as et
						qprint("running with ElementTree")
					except ImportError:
						eprint("Fatal error: Failed to import ElementTree from any known place. XML writing is impossible. ")
	############################################################################################################################################################################################
	# On with bulk of code
	old_time=start_time # for speed diagnostics
	for filename in args.files:
		# if input file doesn't have _pro suffix, rename input files to have this suffix. Later will write data if it is missing.
		processed_suffix='_pro'
		if filename[-7:-3]!='_pro':
			processed_suffix='' #temporary fix
		hf = h5py.File(filename,'r')
		dims=[]
		for dim in hf['mesh']['array_dimensions']: #full dimensions of mesh
			dims.append(dim) 
		extents = []  #IE the dimensions actually used in this particular run
		extents.append(hf['mesh']['radial_index_bound'][1])
		extents.append(hf['mesh']['theta_index_bound'][1])
		extents.append(hf['mesh']['phi_index_bound'][1])
		if args.dimensions:
			extents[0]=args.dimensions[0]	
		entities={
			"Extent_r":str(extents[0]),
			"Slice_YZ_extent_theta":str(extents[1]*2),
			"Slice_XZ_extent_theta":str(extents[2]*2),
			"Slice_XY_extent_theta":str(extents[2]),
			"Dim_r":str(extents[0]+1),
		}
		con_sp={
		'X':'YZ',
		'Y':'XZ',
		'Z':'XY'
		}
		topo_type='2DRectMesh'
		geom_type='VXVY'
		is_3d=True
		is_2d=False
		if extents[2]==1:
			# Hopefully the end case will handle both.
			eprint('Error: Tried to run on a 2d file instead of 3d file')
			break
		n_hyperslabs = hf['mesh']['nz_hyperslabs'].value
		n_elemental_species = hf['abundance']['xn_c'].shape[3]
		# explode filename into list
		file_directory=''
		if re.search('.*\/(?!.+\/)',filename):
			file_directory = re.search('.*\/(?!.+\/)',filename).group()
		filename_part = re.search('(?!.*\/).*',filename).group().rsplit('_')
		if args.prefix:
			filename_part[0]=args.prefix
		if args.dir:
			file_directory = str(args.dir)
			qprint("Directory set to "+file_directory)
		if not file_directory.endswith('/') and file_directory != '':
			file_directory+='/'
		extension='.xmf'
		new_filename=[]
		reduced_hf=[]
		if args.xdmf:
			extension='.xdmf'
		if args.shortfilename:
			short=True
		else:
			short=False
		xdmf_filename='2D_'+filename_part[0]+['_step-','-'][short]+filename_part[1]+extension
		file_out_name=file_directory+xdmf_filename
		new_filename={}
		reduced_hf={}
		coords=['X','Y','Z']
		for axis in coords:
			new_filename[axis]='2D_'+filename_part[0]+['_axis-','-'][short]+con_sp[axis]+['_step-','-'][short]+filename_part[1]+'.h5'
			reduced_hf[axis]=h5py.File('2D_'+filename_part[0]+['_axis-','-'][short]+con_sp[axis]+['_step-','-'][short]+filename_part[1]+'.h5','w')
		del axis,short
		######### Module: Mesh setup ###########
		# Here we write out mesh variables for each slice and make corresponding xdmf
		
		# make variable that the normal parser can easily check for and refer to this script
		for coord in coords:
			reduced_hf[coord].create_dataset('is_reduced',data=True)
 		# X slice mesh:
		reduced_hf['X'].create_group('/mesh')
		reduced_hf['X'].create_dataset('/mesh/x_ef',data=hf['mesh']['x_ef'])
		y=np.append(hf['mesh']['y_ef'],hf['mesh']['y_ef'].value[1:]+np.pi)
		reduced_hf['X'].create_dataset('/mesh/y_ef',data=y)
		# Y slice mesh:
		reduced_hf['Y'].create_group('/mesh')
		reduced_hf['Y'].create_dataset('/mesh/x_ef',data=hf['mesh']['x_ef'])
		reduced_hf['Y'].create_dataset('/mesh/y_ef',data=y)
		# del y
		# Z slice mesh:
		reduced_hf['Z'].create_group('/mesh')
		reduced_hf['Z'].create_dataset('/mesh/x_ef',data=hf['mesh']['x_ef'])
		reduced_hf['Z'].create_dataset('/mesh/z_ef',data=hf['mesh']['z_ef'].value)
		# create corresponding xdmf:
		xdmf = et.Element("Xdmf", Version="2.0")
		domain = et.SubElement(xdmf,"Domain")
		grid = {'YZ/Hydro':et.SubElement(domain,"Grid",Name="YZ/Hydro"),
				'XZ/Hydro':et.SubElement(domain,"Grid",Name="XZ/Hydro"),
				'XY/Hydro':et.SubElement(domain,"Grid",Name="XY/Hydro"), #I would have this later but apparently order matters because VisIt
				'YZ/Abundance':et.SubElement(domain,"Grid",Name="YZ/Abundance"),
				'YZ/Radiation':et.SubElement(domain,"Grid",Name="YZ/Radiation"),
				'XZ/Abundance':et.SubElement(domain,"Grid",Name="XZ/Abundance"),
				'XZ/Radiation':et.SubElement(domain,"Grid",Name="XZ/Radiation"),
				'XY/Abundance':et.SubElement(domain,"Grid",Name="XY/Abundance"),
				'XY/Radiation':et.SubElement(domain,"Grid",Name="XY/Radiation"),
				}
		# YZ/Hydro,XZ/Hydro geometry and topology xdmf:
		for name in ['YZ/Hydro','XZ/Hydro']:
			et.SubElement(grid[name],"Topology",TopologyType=topo_type,NumberOfElements=str(extents[1]*2+1)+' &Dim_r;')
			geometry = et.SubElement(grid[name],"Geometry",GeometryType=geom_type)
			for n,coord_name in enumerate(["x_ef","y_ef"]):
				parent_element=geometry
				if coord_name=='x_ef':
					unit_changing_function = et.SubElement(geometry,"DataItem",Dimensions='&Dim_r;',ItemType="Function",Function="$0/100000")
					parent_element=unit_changing_function
				hyperslab = et.SubElement(parent_element, "DataItem",Dimensions=str(['&Dim_r;',extents[n]*2+1][coord_name=='y_ef']),ItemType="HyperSlab")
				et.SubElement(hyperslab,"DataItem",Dimensions="3 1",Format="XML").text="0 1 "+str(['&Dim_r;',extents[n]*2+1][coord_name=='y_ef'])
				et.SubElement(hyperslab,"DataItem",Dimensions=str(dims[n]*[1,2][coord_name=='y_ef']+1),NumberType="Float",Precision="8",Format="HDF").text = "&h5path"+name[0]+";:/mesh/" + coord_name
			et.SubElement(grid[name],"Time",Value=str(hf['mesh']['time'].value-hf['mesh']['t_bounce'].value))
		# For XY/Hydro geometry and topology xdmf:
		et.SubElement(grid['XY/Hydro'],"Topology",TopologyType=topo_type,NumberOfElements=str(extents[2]+1)+' &Dim_r;')
		geometry = et.SubElement(grid['XY/Hydro'],"Geometry",GeometryType=geom_type)
		for n,coord_name in enumerate(["x_ef","z_ef"]):
			# correct for switch from y_ef to z_ef
			if n==1:
				n=2
			parent_element=geometry
			if coord_name=='x_ef':
				unit_changing_function = et.SubElement(geometry,"DataItem",Dimensions=[str(extents[n]+1),'&Dim_r;'][n==0],ItemType="Function",Function="$0/100000")
				parent_element=unit_changing_function
			hyperslab = et.SubElement(parent_element, "DataItem",Dimensions=[str(extents[n]+1),'&Dim_r;'][n==0],ItemType="HyperSlab")
			et.SubElement(hyperslab,"DataItem",Dimensions="3 1",Format="XML").text="0 1 "+[str(extents[n]+1),'&Dim_r;'][n==0]
			et.SubElement(hyperslab,"DataItem",Dimensions=str(dims[n]+1),NumberType="Float",Precision="8",Format="HDF").text = "&h5pathZ;:/mesh/" + coord_name
		et.SubElement(grid['XY/Hydro'],"Time",Value=str(hf['mesh']['time'].value-hf['mesh']['t_bounce'].value))
		# for rest of grids:
		for name in grid:
			if name[3:]!='Hydro':
				et.SubElement(grid[name],"Topology",Reference="/Xdmf/Domain/Grid[@Name='"+name[0:2]+"/Hydro']/Topology[1]")
				et.SubElement(grid[name],"Geometry",Reference="/Xdmf/Domain/Grid[@Name='"+name[0:2]+"/Hydro']/Geometry[1]")
		######### Module: Z Slicer ###########
		qprint("Slicing through Z axis")
		S_P='Z'
		storage_names = { 
			"Entropy":"entropy",
			"Density":"rho_c",
			"v_radial":"u_c",
			"v_theta":"v_c",
			"v_phi":"w_c",
			"BruntViasala_freq":"wBVMD",
			"Electron_fraction":"ye_c",
			"Gravity_phi":"grav_x_c",
			"Gravity_r":"grav_y_c",
			"Gravity_theta":"grav_z_c",
			"Lepton_fraction":"ylep",
			"Mach_number":"v_csound",
			"Neutrino_Heating_rate":"dudt_nu",
			"Nuclear_Heating_rate":"dudt_nuc",
			"Pressure":"press",
			"Temperature":"t_c",
		}
		reduced_hf[S_P].create_group('/fluid')
		extents_str='&Slice_'+con_sp[S_P]+'_extent_theta; &Extent_r;'
		dims_str=str(dims[2])+" "+str(dims[0])
		for name in storage_names:
			qprint('	Writing '+name+' to dataset in '+S_P+' reduced hdf5 file')
			scalar_quantity=hf['fluid'][storage_names[name]][:,dims[1]/2,:]
			for sl in range(2,n_hyperslabs+1):
				i=sl
				if args.repeat:
					i=1
				temp_hf= h5py.File(re.sub("\d\d\.h5",str(format(i, '02d'))+'.h5',re.sub("\d\d_pro\.h5",str(format(i, '02d'))+'_pro.h5',filename)),'r')
				scalar_quantity=np.vstack((scalar_quantity,temp_hf['fluid'][storage_names[name]][:,dims[1]/2,:]))
			reduced_hf[S_P].create_dataset('/fluid/'+storage_names[name],data=scalar_quantity)
			at = et.SubElement(grid[con_sp[S_P]+'/Hydro'],"Attribute",Name=name,AttributeType="Scalar",Center="Cell",Dimensions=extents_str)
			hyperslab = et.SubElement(at, "DataItem",Dimensions=extents_str,ItemType="HyperSlab")
			et.SubElement(hyperslab,"DataItem",Dimensions="3 2",Format="XML").text="0 0 1 1 "+extents_str
			et.SubElement(hyperslab,"DataItem",Dimensions=dims_str,NumberType="Float",Precision="8",Format="HDF").text = "&h5path"+S_P+";:/fluid/" + storage_names[name]
		#	For abundance:
		at = et.SubElement(grid[con_sp[S_P]+'/Abundance'],"Attribute",Name='nse_flag',AttributeType="Scalar",Center="Cell",Dimensions=extents_str)
		hyperslab = et.SubElement(at, "DataItem",Dimensions=extents_str,ItemType="HyperSlab")
		et.SubElement(hyperslab,"DataItem",Dimensions="3 2",Format="XML").text="0 0 1 1 "+extents_str
		scalar_quantity=hf['abundance']['nse_c'][:,dims[1]/2,:]
		for sl in range(2,n_hyperslabs+1):
			i=sl
			if args.repeat:
				i=1
			temp_hf= h5py.File(re.sub("\d\d\.h5",str(format(i, '02d'))+'.h5',re.sub("\d\d_pro\.h5",str(format(i, '02d'))+'_pro.h5',filename)),'r')
			scalar_quantity=np.vstack((scalar_quantity,temp_hf['abundance']['nse_c'][:,dims[1]/2,:]))
		reduced_hf[S_P].create_dataset('/abundance/nse_flag',data=scalar_quantity)
		et.SubElement(hyperslab,"DataItem",Dimensions=re.sub(r'[^\w\ ]','',str(scalar_quantity.shape)),NumberType="Int",Format="HDF").text = "&h5path"+S_P+";:/abundance/nse_flag"

		scalar_quantity=hf['abundance']['xn_c'][:,dims[1]/2,:,:]
		for sl in range(2,n_hyperslabs+1):
			i=sl
			if args.repeat:
				i=1
			temp_hf= h5py.File(re.sub("\d\d\.h5",str(format(i, '02d'))+'.h5',re.sub("\d\d_pro\.h5",str(format(i, '02d'))+'_pro.h5',filename)),'r')
			scalar_quantity=np.vstack((scalar_quantity,temp_hf['abundance']['xn_c'][:,dims[1]/2,:]))
		reduced_hf[S_P].create_dataset('/abundance/xn_c',data=scalar_quantity)
		species_names=hf['abundance']['a_name'].value
		if n_elemental_species-1==species_names.shape[0]:
			species_names=np.append(species_names,'aux')
		if not args.disable or "abundance" not in args.disable:
			for el,name in enumerate(species_names):
				if re.findall('\D\d',name): #if there is a transition between a non digit to a digit in the element name (IE in "li3" it would match because of the "i3")
					element_name=re.sub('\d','',name).capitalize() #set element_name to the capitalized element without the number
					name=element_name+re.sub('\D','',name) #find the transition between elements name and number
					if not grid.has_key(S_P+'Abundance'+'/'+element_name): #If the grid for that element doesn't already exist, create it 
						grid[con_sp[S_P]+'/Abundance'+'/'+element_name]=et.SubElement(domain,"Grid",Name=con_sp[S_P]+'/Abundance'+'/'+element_name,GridType="Uniform")
						et.SubElement(grid[con_sp[S_P]+'/Abundance'+'/'+element_name],"Topology",Reference="/Xdmf/Domain/Grid[@Name='"+con_sp[S_P]+"/Abundance']/Topology[1]")
						et.SubElement(grid[con_sp[S_P]+'/Abundance'+'/'+element_name],"Geometry",Reference="/Xdmf/Domain/Grid[@Name='"+con_sp[S_P]+"/Abundance']/Geometry[1]")
					attribute=et.SubElement(grid[con_sp[S_P]+'/Abundance'+'/'+element_name],"Attribute",Name=name,AttributeType="Scalar",Center="Cell")
				else:
					attribute=et.SubElement(grid[con_sp[S_P]+'/Abundance'],"Attribute",Name=name,AttributeType="Scalar",Center="Cell")
				dataElement = et.SubElement(attribute,"DataItem", ItemType="HyperSlab", Dimensions=extents_str)
				et.SubElement(dataElement,"DataItem",Dimensions="3 3",Format="XML").text="0 0 "+str(el)+" 1 1 1 "+extents_str+" 1"
				et.SubElement(dataElement,"DataItem",Dimensions=dims_str+" "+str(n_elemental_species),NumberType="Float",Precision="8",Format="HDF").text= "&h5path"+S_P+";:/abundance/xn_c"
				n+=1
		######### Module: X Slicer ###########
		qprint("Slicing through X axis")
		S_P='X' #S_P stands for Slice Prefix
		reduced_hf[S_P].create_group('/fluid')
		extents_str="&Slice_"+con_sp[S_P]+"_extent_theta; &Extent_r;"
		dims_str=str(dims[1]*2)+" "+str(dims[0])
		i=n_hyperslabs/2
		if args.repeat:
			i=1
		hf2=h5py.File(re.sub("\d\d\.h5",str(format(i, '02d'))+'.h5',filename),'r')
		#	For hydro:
		for name in storage_names:
			qprint('	Writing '+name+' to '+con_sp[S_P]+' dataset in reduced hdf5 file')
			scalar_quantity=np.vstack((hf['fluid'][storage_names[name]][0,:,:],hf2['fluid'][storage_names[name]].value[0,::-1,:]))
			reduced_hf[S_P].create_dataset('/fluid/'+storage_names[name],data=scalar_quantity)
			at = et.SubElement(grid[con_sp[S_P]+'/Hydro'],"Attribute",Name=name,AttributeType="Scalar",Center="Cell",Dimensions=extents_str)
			hyperslab = et.SubElement(at, "DataItem",Dimensions=extents_str,ItemType="HyperSlab")
			et.SubElement(hyperslab,"DataItem",Dimensions="3 2",Format="XML").text="0 0 1 1 "+extents_str
			et.SubElement(hyperslab,"DataItem",Dimensions=dims_str,NumberType="Float",Precision="8",Format="HDF").text = "&h5path"+S_P+";:/fluid/" + storage_names[name]
		#	For abundance:

		at = et.SubElement(grid[con_sp[S_P]+'/Abundance'],"Attribute",Name='nse_flag',AttributeType="Scalar",Center="Cell",Dimensions=extents_str)
		hyperslab = et.SubElement(at, "DataItem",Dimensions=extents_str,ItemType="HyperSlab")
		et.SubElement(hyperslab,"DataItem",Dimensions="3 2",Format="XML").text="0 0 1 1 "+extents_str
		scalar_quantity=np.vstack((hf['abundance']['nse_c'][0,:,:],hf2['abundance']['nse_c'].value[0,::-1,:]))
		reduced_hf[S_P].create_dataset('/abundance/nse_flag',data=scalar_quantity)
		et.SubElement(hyperslab,"DataItem",Dimensions=re.sub(r'[^\w\ ]','',str(scalar_quantity.shape)),NumberType="Int",Format="HDF").text = "&h5path"+S_P+";:/abundance/nse_flag"
		
		scalar_quantity=np.vstack((hf['abundance']['xn_c'][0,:,:,:],hf2['abundance']['xn_c'].value[0,::-1,:,:]))
		reduced_hf[S_P].create_dataset('/abundance/xn_c',data=scalar_quantity)
		species_names=hf['abundance']['a_name'].value
		if n_elemental_species-1==species_names.shape[0]:
			species_names=np.append(species_names,'aux')
		if not args.disable or "abundance" not in args.disable:
			for el,name in enumerate(species_names):
				if re.findall('\D\d',name): #if there is a transition between a non digit to a digit in the element name (IE in "li3" it would match because of the "i3")
					element_name=re.sub('\d','',name).capitalize() #set element_name to the capitalized element without the number
					name=element_name+re.sub('\D','',name) #find the transition between elements name and number
					if not grid.has_key(con_sp[S_P]+'/Abundance'+'/'+element_name): #If the grid for that element doesn't already exist, create it 
						grid[con_sp[S_P]+'/Abundance/'+element_name]=et.SubElement(domain,"Grid",Name=con_sp[S_P]+'/Abundance/'+element_name,GridType="Uniform")
						et.SubElement(grid[con_sp[S_P]+'/Abundance/'+element_name],"Topology",Reference="/Xdmf/Domain/Grid[@Name='"+con_sp[S_P]+"/Abundance']/Topology[1]")
						et.SubElement(grid[con_sp[S_P]+'/Abundance/'+element_name],"Geometry",Reference="/Xdmf/Domain/Grid[@Name='"+con_sp[S_P]+"/Abundance']/Geometry[1]")
					attribute=et.SubElement(grid[con_sp[S_P]+'/Abundance'+'/'+element_name],"Attribute",Name=name,AttributeType="Scalar",Center="Cell")
				else:
					attribute=et.SubElement(grid[con_sp[S_P]+'/Abundance'],"Attribute",Name=name,AttributeType="Scalar",Center="Cell")
				dataElement = et.SubElement(attribute,"DataItem", ItemType="HyperSlab", Dimensions=extents_str)
				et.SubElement(dataElement,"DataItem",Dimensions="3 3",Format="XML").text="0 0 "+str(el)+" 1 1 1 "+extents_str+" 1"
				et.SubElement(dataElement,"DataItem",Dimensions=dims_str+" "+str(n_elemental_species),NumberType="Float",Precision="8",Format="HDF").text= "&h5path"+S_P+";:/abundance/xn_c"
				n+=1
		hf2.close()
		######### Module: Y Slicer ###########
		qprint("Slicing through Y axis")
		S_P='Y'
		i=n_hyperslabs/4
		if args.repeat:
			i=1
		hf2=h5py.File(re.sub("\d\d\.h5",str(format(i, '02d'))+'.h5',filename),'r')
		reduced_hf[S_P].create_group('/fluid')
		i=n_hyperslabs*3/4
		if args.repeat:
			i=1
		hf3=h5py.File(re.sub("\d\d\.h5",str(format(i, '02d'))+'.h5',filename),'r')
		# for Hydro:

		for name in storage_names:
			qprint('	Writing '+name+' to '+con_sp[S_P]+' dataset in reduced hdf5 file')
			scalar_quantity=np.vstack((hf2['fluid'][storage_names[name]][0,:,:],hf3['fluid'][storage_names[name]].value[0,::-1,:]))
			reduced_hf[S_P].create_dataset('/fluid/'+storage_names[name],data=scalar_quantity)
			at = et.SubElement(grid[con_sp[S_P]+'/Hydro'],"Attribute",Name=name,AttributeType="Scalar",Center="Cell",Dimensions=extents_str)
			hyperslab = et.SubElement(at, "DataItem",Dimensions=extents_str,ItemType="HyperSlab")
			et.SubElement(hyperslab,"DataItem",Dimensions="3 2",Format="XML").text="0 0 1 1 "+extents_str
			et.SubElement(hyperslab,"DataItem",Dimensions=dims_str,NumberType="Float",Precision="8",Format="HDF").text = "&h5path"+S_P+";:/fluid/" + storage_names[name]
		#	for Abundance:

		at = et.SubElement(grid[con_sp[S_P]+'/Abundance'],"Attribute",Name='nse_flag',AttributeType="Scalar",Center="Cell",Dimensions=extents_str)
		hyperslab = et.SubElement(at, "DataItem",Dimensions=extents_str,ItemType="HyperSlab")
		et.SubElement(hyperslab,"DataItem",Dimensions="3 2",Format="XML").text="0 0 1 1 "+extents_str
		scalar_quantity=np.vstack((hf2['abundance']['nse_c'][0,:,:],hf3['abundance']['nse_c'].value[0,::-1,:]))
		reduced_hf[S_P].create_dataset('/abundance/nse_flag',data=scalar_quantity)
		et.SubElement(hyperslab,"DataItem",Dimensions=re.sub(r'[^\w\ ]','',str(scalar_quantity.shape)),NumberType="Int",Format="HDF").text = "&h5path"+S_P+";:/abundance/nse_flag"
		

		scalar_quantity=np.vstack((hf2['abundance']['xn_c'][0,:,:,:],hf3['abundance']['xn_c'].value[0,::-1,:,:]))
		reduced_hf[S_P].create_dataset('/abundance/xn_c',data=scalar_quantity)
		species_names=hf['abundance']['a_name'].value
		if n_elemental_species-1==species_names.shape[0]:
			species_names=np.append(species_names,'aux')
		if not args.disable or "abundance" not in args.disable:
			for el,name in enumerate(species_names):
				if re.findall('\D\d',name): #if there is a transition between a non digit to a digit in the element name (IE in "li3" it would match because of the "i3")
					element_name=re.sub('\d','',name).capitalize() #set element_name to the capitalized element without the number
					name=element_name+re.sub('\D','',name) #find the transition between elements name and number
					if not grid.has_key(con_sp[S_P]+'/Abundance/'+element_name): #If the grid for that element doesn't already exist, create it 
						grid[con_sp[S_P]+'/Abundance'+'/'+element_name]=et.SubElement(domain,"Grid",Name=con_sp[S_P]+'/Abundance/'+element_name,GridType="Uniform")
						et.SubElement(grid[con_sp[S_P]+'/Abundance'+'/'+element_name],"Topology",Reference="/Xdmf/Domain/Grid[@Name='"+con_sp[S_P]+"/Abundance']/Topology[1]")
						et.SubElement(grid[con_sp[S_P]+'/Abundance'+'/'+element_name],"Geometry",Reference="/Xdmf/Domain/Grid[@Name='"+con_sp[S_P]+"/Abundance']/Geometry[1]")
					attribute=et.SubElement(grid[con_sp[S_P]+'/Abundance'+'/'+element_name],"Attribute",Name=name,AttributeType="Scalar",Center="Cell")
				else:
					attribute=et.SubElement(grid[con_sp[S_P]+'/Abundance'],"Attribute",Name=name,AttributeType="Scalar",Center="Cell")
				dataElement = et.SubElement(attribute,"DataItem", ItemType="HyperSlab", Dimensions=extents_str)
				et.SubElement(dataElement,"DataItem",Dimensions="3 3",Format="XML").text="0 0 "+str(el)+" 1 1 1 "+extents_str+" 1"
				et.SubElement(dataElement,"DataItem",Dimensions=dims_str+" "+str(n_elemental_species),NumberType="Float",Precision="8",Format="HDF").text= "&h5path"+S_P+";:/abundance/xn_c"
				n+=1
		hf2.close()
		hf3.close()
		######### Module: Auxilary Variables ##############################################################################################################################################################
		if args.aux:
			qprint("Creating derived values")
			n_groups=hf['radiation']['raddim'][0]
			n_species=hf['radiation']['raddim'][1]
			energy_edge=hf['radiation']['unubi'].value
			energy_center=hf['radiation']['unui'].value
			d_energy=[]
			for i in range(0,n_groups):
				d_energy.append(energy_edge[i+1]-energy_edge[i])
			d_energy=np.array(d_energy)
			e3de = energy_center**3*d_energy
			e5de = energy_center**5*d_energy
			cvel=2.99792458e10
			pi=3.141592653589793
			ergmev = 1.602177e-6
			h = 4.13567e-21
			ecoef = 4.0 * pi * ergmev/(h*cvel)**3
			######## Compute E_RMS_array (size N_species) of arrays (size N_groups) ##############
			# # initialize variables for parallel loop
			num_cores=min(16,mp.cpu_count())
			E_RMS_array=np.empty((n_hyperslabs,n_species,dims[2],dims[0]))
			def compute_E_RMS_array_z(sl):	
				sl+=1
				i=sl
				if args.repeat:
					i=1
				qprint("	Computing E_RMS_[1.."+str(n_species)+"] for slice "+str(sl)+" from "+re.sub("\d\d\.h5",str(format(i, '02d'))+'.h5',re.sub("\d\d_pro\.h5",str(format(i, '02d'))+'_pro.h5',filename)))
				temp_hf= h5py.File(re.sub("\d\d\.h5",str(format(i, '02d'))+'.h5',re.sub("\d\d_pro\.h5",str(format(i, '02d'))+'_pro.h5',filename)),'r')
				psi0_c=temp_hf['radiation']['psi0_c'][:,dims[1]/2,:,:]
				row=np.empty((n_species,dims[2]/n_hyperslabs,dims[0]))
				for n in range(0,n_species):
					numerator=np.sum(psi0_c[:,:,n]*e5de,axis=2)
					denominator=np.sum(psi0_c[:,:,n]*e3de,axis=2)
					row[n][:][:]=np.sqrt(numerator/(denominator+1e-100))
				return row
			def compute_E_RMS_array_x():	
				i=n_hyperslabs
				if args.repeat:
					i=1
				qprint("	Computing E_RMS_[1.."+str(n_species)+"] for slice "+str(sl)+" from "+re.sub("\d\d\.h5",str(format(i, '02d'))+'.h5',re.sub("\d\d_pro\.h5",str(format(i, '02d'))+'_pro.h5',filename)))
				temp_hf= h5py.File(re.sub("\d\d\.h5",str(format(i, '02d'))+'.h5',re.sub("\d\d_pro\.h5",str(format(i, '02d'))+'_pro.h5',filename)),'r')
				psi0_c=hf['radiation']['psi0_c'][0,:,:,:]
				row=np.empty((n_species,dims[2]*2,dims[0]))
				for n in range(0,n_species):
					numerator=np.sum(psi0_c[:,:,n,:]*e5de,axis=2)
					denominator=np.sum(psi0_c[:,:,n,:]*e3de,axis=2)
					row[n][0:dims[2]][:]=np.sqrt(numerator/(denominator+1e-100))
				psi0_c=hf['radiation']['psi0_c'][0,:,:,:]
				for n in range(0,n_species):
					numerator=np.sum(psi0_c[:,:,n,:]*e5de,axis=2)
					denominator=np.sum(psi0_c[:,:,n,:]*e3de,axis=2)
					row[n][dims[2]:dims[2]*2][:]=np.sqrt(numerator/(denominator+1e-100))
				return row.shape
			def radiation_xdmf(xdmf_alias_name,hdf_variable_name,S_P,xdmf_grid,hdf_directory,dims_str,extents_str):
				# Fix potential directory issue
				if hdf_directory[0]!='/':
					hdf_directory='/'+hdf_directory
				if hdf_directory[:-1]!='/':
					hdf_directory+='/'
				if xdmf_grid[0]!='/':
					xdmf_grid='/'+xdmf_grid
				# Make xdmf elements
				at = et.SubElement(grid[con_sp[S_P]+xdmf_grid],"Attribute",Name=xdmf_alias_name,AttributeType="Scalar",Center="Cell",Dimensions=extents_str)
				hyperslab = et.SubElement(at, "DataItem",Dimensions=extents_str,ItemType="HyperSlab")
				et.SubElement(hyperslab,"DataItem",Dimensions="3 2",Format="XML").text="0 0 1 1 "+extents_str
				et.SubElement(hyperslab,"DataItem",Dimensions=dims_str,NumberType="Float",Precision="8",Format="HDF").text = "&h5path"+S_P+";:"+hdf_directory+hdf_variable_name
			S_P='Z'
			qprint("Z slice E_RMS:")
			results = np.array((n_hyperslabs,n_species,dims[2],dims[0]))
			results = Parallel(n_jobs=num_cores)(delayed(compute_E_RMS_array_z)(sl) for sl in range(0,n_hyperslabs))
			#concatenate together the member of each array within E_RMS_ARRAY and write to auxilary HDF file
			qprint("Concatenating E_RMS_[0.."+str(n_species)+"] results...")
			E_RMS_array=np.hstack(results)
			qprint("Done.\nWriting results:")
			for n,sp in enumerate(['e','e-bar','mt','mt-bar']):
				qprint("	Writing /Radiation/E_RMS_"+sp+" out to "+con_sp[S_P]+" HDF5 file")
				reduced_hf[S_P].create_dataset("/radiation/E_RMS_"+sp,data=E_RMS_array[n])
			del E_RMS_array, results
			# Corresponding xdmf:
			extents_str="&Slice_"+con_sp[S_P]+"_extent_theta; &Extent_r;"
			dims_str=str(dims[2])+" "+str(dims[0])
			for n in ['e','e-bar','mt','mt-bar']:
				qprint("	Creating "+con_sp[S_P]+"/Radiation/E_RMS_"+str(n)+" xdmf element")
				radiation_xdmf('E_RMS_'+str(n),'E_RMS_'+str(n),S_P,'/Radiation','/radiation',dims_str,extents_str)

			extents_str=str(extents[1]*2)+" "+str(extents[0])
			dims_str=str(dims[1]*2)+" "+str(dims[0])
			# ######## luminosity part ###########
			# qprint("Computing luminosities")
			# psi1_e=hf['radiation']['psi1_e']
			# radius=hf['mesh']['x_ef'].value
			# agr_e=hf['fluid']['agr_e'].value
			# cell_area_GRcorrected=4*pi*radius**2/agr_e**4
			# lumin=np.empty((dims[2],dims[1],dims[0]))
			# lumin_array=[]
			# def compute_luminosity(species,sl):
			# 	j=sl
			# 	if args.repeat:
			# 		sl=0
			# 	temp_hf= h5py.File(re.sub("\d\d\.h5",str(format(sl+1, '02d'))+'.h5',re.sub("\d\d_pro\.h5",str(format(sl+1, '02d'))+'_pro.h5',filename)),'r')
			# 	psi1_e=temp_hf['radiation']['psi1_e']
			# 	qprint("	On slice "+str(j+1)+" of "+str(n_hyperslabs)+" from "+re.sub("\d\d\.h5",str(format(sl+1, '02d'))+'.h5',re.sub("\d\d_pro\.h5",str(format(sl+1, '02d'))+'_pro.h5',filename)))
			# 	br()
			# 	del temp_hf
			# 	return np.sum(psi1_e[:,:,:,species]*e3de, axis=3)*np.tile(cell_area_GRcorrected[1:dims[0]+1],(dims[2]/n_hyperslabs,dims[1],1))*(cvel*ecoef*1e-51)
			# lumin_array=np.empty((n_species,dims[2],dims[1],dims[0]))
			# for species in range(0,n_species):
			# 	qprint("Computing luminosity for species "+str(species)+":")
			# 	lumin_array[species] = np.vstack(Parallel(n_jobs=num_cores)(delayed(compute_luminosity)(0,sl) for sl in range(0,n_hyperslabs)))
			# qprint("########################################")
			# for n,lumin in enumerate(lumin_array):
			# 	qprint("Writing luminosity species "+str(n)+" to auxilary hdf file")
			# 	reduced_hf.create_dataset("/radiation/Luminosity_"+str(n),data=lumin)
			# del lumin_array,lumin,mask
		for S_P in coords:
			reduced_hf[S_P].close()
		############################################################################################################################################################################################
		######### Module: Xdmf Document writer ###########
		try:
			f=open(file_out_name,'w')
			del extension,file_directory
			# if lxml module loaded use it to write document (fasted, simplest implementation):
			for coord in coords:
				entities['h5path'+coord]=new_filename[coord]
			entity_str = ''
			for key,value in six.iteritems(entities):
				entity_str+="\n  <!ENTITY "+key+" \""+value+"\">"
			entity_str+="\n  <!-- Note that Dim_r must be exactly 1 more than Extent_r or VisIt will have a spontanious freak out session -->"
			try:
				#write to file:
				f.write(\
					#remove all the '&amp;' tags that appear due to xml parser and replace with '&' so the aliasing works
					re.sub(\
						'&amp;','&',et.tostring(xdmf,\
							pretty_print=True,\
							xml_declaration=True,\
							encoding="ASCII",\
							doctype="<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" ["+entity_str+"\n]>"\
							)\
					)\
				)
			#other ElementTree writers can use this slower writer that does the same thing:
			except:
				f.close()
				import xml.etree.cElementTree as et
				import xml.dom.minidom as md
				f=open(file_out_name,'w')
				qprint("Writing "+file_out_name+" with improvised \"pretty print\"")
				def prettify(elem):
					rough_string = et.tostring(elem, 'ASCII')
					reparsed = md.parseString(rough_string)
					t = re.sub('&amp;','&',reparsed.toprettyxml(indent="  "))#lxml prettyprint uses 2 spaces per step so we do here as well
					t = ''.join(t.splitlines(True)[1:]) #removes extra doc declaration that mysteriously appears
					return t
				# write custom doctype declaration
				f.write("<?xml version='1.0' encoding='ASCII'?>\n<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" ["+entity_str+"\n]>\n")
				f.close()
				f=open(file_out_name,'a')
				f.write(prettify(xdmf))
			f.close()
			qprint("--- "+file_out_name+" created in %s seconds ---" % (time.time()-old_time))
			old_time=time.time()
		except:
			eprint("Fatal error:")
			eprint("	File write error. Try adding adding a br() (local name of pdb.set_trace()) in the writer section at the end of script to debug")
			sys.exit()