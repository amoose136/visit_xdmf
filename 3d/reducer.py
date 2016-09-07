#!/usr/bin/python
#  ______   ______  _____    _    _   ______  ______  ______  
# | |  | \ | |     | | \ \  | |  | | | |     | |     | |  | \ 
# | |__| | | |---- | |  | | | |  | | | |     | |---- | |__| | 
# |_|  \_\ |_|____ |_|_/_/  \_|__|_| |_|____ |_|____ |_|  \_\ 
# 		A script to reduce 3d datasets to 2d ones                
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
	# parser.add_argument('--extents','-e',dest='dimensions',metavar='int',action='store',type=int,nargs=3, help='dimensions to crop (by grid cell number not spatial dimensions)')
	parser.add_argument('--slices',dest='slices',metavar='int',action='store',type=int,nargs='?', help='number of slices to use')
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
		topo_type='2DRectMesh'
		geom_type='VXVY'
		is_3d=True
		is_2d=False
		if extents[2]==1:
			eprint('Error: Tried to run on a 2d file instead of 3d file')
			break
		n_hyperslabs = hf['mesh']['nz_hyperslabs'].value
		n_elemental_species = hf['abundance']['xn_c'].shape[3]
		if args.slices:
			if not args.slices>slices:
				slices=args.slices
			else:
				eprint("Error: slices must not be more than the number of wedges")
				sys.exit()
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
		if args.xdmf:
			extension='.xdmf'
		if args.shortfilename:
			XDMF_filename='2D_'+filename_part[0]+'-'+filename_part[3]+'-'+filename_part[1]+extension
		else:
			XDMF_filename='2D_'+filename_part[0]+'_grid-'+filename_part[3]+'_step-'+filename_part[1]+extension
		file_out_name=file_directory+XDMF_filename
		new_filename=XDMF_filename[:-len(extension)]+'.h5'
		reduced_hf=h5py.File(XDMF_filename[:-len(extension)]+'.h5','w')
		######### Module: Mesh setup ###########
		# Here we write out mesh variables for each slice and make corresponding xdmf
		
		# make variable that the normal parser can easily check for and refer to this script
		reduced_hf.create_dataset('/is_reduced',data=True)
		# X slice mesh:
		reduced_hf.create_group('/X')
		reduced_hf.create_group('/X/mesh')
		reduced_hf.create_dataset('/X/mesh/x_ef',data=hf['mesh']['x_ef'])
		y=np.append(hf['mesh']['y_ef'],hf['mesh']['y_ef'].value[1:]+np.pi)
		reduced_hf.create_dataset('/X/mesh/y_ef',data=y)
		# Y slice mesh:
		reduced_hf.create_group('/Y')
		reduced_hf.create_group('/Y/mesh')
		reduced_hf.create_dataset('/Y/mesh/x_ef',data=hf['mesh']['x_ef'])
		reduced_hf.create_dataset('/Y/mesh/y_ef',data=y)
		# del y
		# Z slice mesh:
		reduced_hf.create_group('/Z')
		reduced_hf.create_group('/Z/mesh')
		reduced_hf.create_dataset('/Z/mesh/x_ef',data=hf['mesh']['x_ef'])
		reduced_hf.create_dataset('/Z/mesh/z_ef',data=hf['mesh']['z_ef'])
		# create corresponding xdmf:
		xdmf = et.Element("Xdmf", Version="2.0")
		domain = et.SubElement(xdmf,"Domain")
		grid = {'X/Hydro':et.SubElement(domain,"Grid",Name="X/Hydro"),
				'X/Abundance':et.SubElement(domain,"Grid",Name="X/Abundance"),
				'X/Radiation':et.SubElement(domain,"Grid",Name="X/Radiation"),
				'X/Mesh':et.SubElement(domain,"Grid",Name="X/Mesh"),
				'Y/Hydro':et.SubElement(domain,"Grid",Name="Y/Hydro"),
				'Y/Abundance':et.SubElement(domain,"Grid",Name="Y/Abundance"),
				'Y/Radiation':et.SubElement(domain,"Grid",Name="Y/Radiation"),
				'Y/Mesh':et.SubElement(domain,"Grid",Name="Y/Mesh"),
				'Z/Hydro':et.SubElement(domain,"Grid",Name="Z/Hydro"),
				'Z/Abundance':et.SubElement(domain,"Grid",Name="Z/Abundance"),
				'Z/Radiation':et.SubElement(domain,"Grid",Name="Z/Radiation"),
				'Z/Mesh':et.SubElement(domain,"Grid",Name="Z/Mesh")
				}
		# X/Hydro geometry and topology:
		et.SubElement(grid['X/Hydro'],"Topology",TopologyType=topo_type,NumberOfElements=str(extents[1]*2+1)+' '+str(extents[0]+1))
		geometry = et.SubElement(grid['X/Hydro'],"Geometry",GeometryType=geom_type)
		coords=["x_ef","y_ef"]
		for n,coord_name in enumerate(coords):
			parent_element=geometry
			if coord_name=='x_ef':
				unit_changing_function = et.SubElement(geometry,"DataItem",Dimensions=str(extents[n]+1),ItemType="Function",Function="$0/100000")
				parent_element=unit_changing_function
			hyperslab = et.SubElement(parent_element, "DataItem",Dimensions=str(extents[n]*[1,2][coord_name=='y_ef']+1),ItemType="HyperSlab")
			et.SubElement(hyperslab,"DataItem",Dimensions="3 1",Format="XML").text="0 1 "+str(extents[n]*[1,2][coord_name=='y_ef']+1)
			et.SubElement(hyperslab,"DataItem",Dimensions=str(dims[n]*[1,2][coord_name=='y_ef']+1),NumberType="Float",Precision="8",Format="HDF").text = "&h5path;:X/mesh/" + coord_name
		et.SubElement(grid['X/Hydro'],"Time",Value=str(hf['mesh']['time'].value-hf['mesh']['t_bounce'].value))
		# For the rest of the X meshes and Y meshes:
		for name in grid:
			if name!='X/Hydro' and name[0]!='Z':
				et.SubElement(grid[name],"Topology",Reference="/Xdmf/Domain/Grid[@Name='X/Hydro']/Topology[1]")
				et.SubElement(grid[name],"Geometry",Reference="/Xdmf/Domain/Grid[@Name='X/Hydro']/Geometry[1]")
		# For Z/Hydro mesh:
		et.SubElement(grid['Z/Hydro'],"Topology",TopologyType=topo_type,NumberOfElements=str(extents[2]+1)+' '+str(extents[0]+1))
		geometry = et.SubElement(grid['Z/Hydro'],"Geometry",GeometryType=geom_type)
		coords=["x_ef","z_ef"]
		for n,coord_name in enumerate(coords):
			# correct for switch from y_ef to z_ef
			if n==1:
				n=2
			parent_element=geometry
			if coord_name=='x_ef':
				unit_changing_function = et.SubElement(geometry,"DataItem",Dimensions=str(extents[n]+1),ItemType="Function",Function="$0/100000")
				parent_element=unit_changing_function
			hyperslab = et.SubElement(parent_element, "DataItem",Dimensions=str(extents[n]+1),ItemType="HyperSlab")
			et.SubElement(hyperslab,"DataItem",Dimensions="3 1",Format="XML").text="0 1 "+str(extents[n]+1)
			et.SubElement(hyperslab,"DataItem",Dimensions=str(dims[n]+1),NumberType="Float",Precision="8",Format="HDF").text = "&h5path;:Z/mesh/" + coord_name
		# for rest of Z grids:
		for name in grid:
			if name!='Z/Hydro' and name[0]=='Z':
				et.SubElement(grid[name],"Topology",Reference="/Xdmf/Domain/Grid[@Name='Z/Hydro']/Topology[1]")
				et.SubElement(grid[name],"Geometry",Reference="/Xdmf/Domain/Grid[@Name='Z/Hydro']/Geometry[1]")
		######### Module: Z Slicer ###########
		# Storage_names dictionary defines mapping between names for scalars as they appear in VisIt (key) and how they are defined in the h5 file (value)
		storage_names = { 
			"Entropy":"entropy",
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
			# "mask":"",#computed quantity to be added later
			"nse_flag":"e_int",
		}
		reduced_hf.create_group('/Z/fluid')
		extents_str=str(extents[1])+" "+str(extents[0])
		dims_str=str(dims[1])+" "+str(dims[0])
		for name in storage_names:
			qprint('Writing '+name+' to reduced dataset hdf5 file')
			scalar_quantity=hf['fluid'][storage_names[name]][:,dims[1]/2,:]
			for sl in range(2,n_hyperslabs+1):
				i=sl
				if args.repeat:
					i=1
				temp_hf= h5py.File(re.sub("\d\d\.h5",str(format(i, '02d'))+'.h5',re.sub("\d\d_pro\.h5",str(format(i, '02d'))+'_pro.h5',filename)),'r')
				scalar_quantity=np.vstack((scalar_quantity,temp_hf['fluid'][storage_names[name]][:,dims[1]/2,:]))
			reduced_hf.create_dataset('/Z/fluid/'+storage_names[name],data=scalar_quantity)
			at = et.SubElement(grid['Z/Hydro'],"Attribute",Name=name,AttributeType="Scalar",Center="Cell",Dimensions=extents_str)
			hyperslab = et.SubElement(at, "DataItem",Dimensions=extents_str,ItemType="HyperSlab")
			et.SubElement(hyperslab,"DataItem",Dimensions="3 2",Format="XML").text="0 0 1 1 "+extents_str
			et.SubElement(hyperslab,"DataItem",Dimensions=dims_str,NumberType="Float",Precision="8",Format="HDF").text = "&h5path;:Z/fluid/" + storage_names[name]

		reduced_hf.create_group('/Z/abundance')
		# Write document tree to file
		try:
			f=open(file_out_name,'w')
			del extension,file_directory
			# if lxml module loaded use it to write document (fasted, simplest implementation):
			entities={
			"h5path":new_filename
			}

			entity_str = ''
			for entity in entities:
				entity_str+="\n  "+entity
			try:
				#write to file:
				f.write(\
					#remove all the '&amp;' tags that appear due to xml parser and replace with '&' so the aliasing works
					re.sub(\
						'&amp;','&',et.tostring(xdmf,\
							pretty_print=True,\
							xml_declaration=True,\
							encoding="ASCII",\
							doctype="<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" [\n  <!ENTITY h5path \""+entities['h5path']+"\">\n]>"\
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
				f.write("<?xml version='1.0' encoding='ASCII'?>\n<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" [\n  <!ENTITY h5path \""+entities['h5path']+"\">\n]>\n")
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



