#!/usr/bin/env python
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
	# sys.path.append('/lustre/atlas/proj-shared/ast109/amos/lib/python2.7/site-packages')
	sys.path.append('./lib/python2.7/site-packages') #added for the assumption that script is called from base directory
	sys.path.append('/sw/redhat6/visit/current/linux-x86_64/lib/site-packages/')
if socket.gethostname()[:5]=='titan':
	sys.path.append('/lustre/atlas/proj-shared/ast109/amos/lib/python2.7/site-packages')
# For my laptop
elif socket.gethostname()=='Lycoris':
	sys.path.append('/Applications/VisIt.app/Contents/Resources/2.10.2/darwin-x86_64/lib/site-packages')

import six
if __name__ == '__main__':
	# try:
	# 	import visit
	# 	import visit_utils
	# except ImportError:
	# 	eprint("Error: visit module import failed\n")
	# 	eprint('Please make sure visit\'s path is in $PATH')
	# 	eprint('	( /path/to/visit/bin )')
	# 	eprint('Please make sure visit\'s /lib/site-packages directory is in $PYTHONPATH')
	# 	eprint ('	( /path/to/visit/VersionNumber/platform/lib/site-packages )')
	# 	sys.exit()
	############################################################################################################################################################################################
	# Contstruct Parser:
	parser = argparse.ArgumentParser(description="Generate XDMF files from Chimera hdf5 files")
	# parser.files is a list of 1 or more h5 files that will have xdmf files generated for them
	parser.add_argument('files',metavar='foo.h5',type=str,nargs='+',help='hdf5 files to process (1 or more args)')
	# parser.add_argument('--extents','-e',dest='dimensions',metavar='int',action='store',type=int,nargs=3, help='dimensions to crop (by grid cell number not spatial dimensions)')
	parser.add_argument('--threads',dest='threads',action='store',metavar='int',type=int, nargs=1, help='specify number of threads')
	parser.add_argument('--slices',dest='slices',metavar='int',action='store',type=int,nargs='?', help='number of slices to use')
	parser.add_argument('--prefix','-p',dest='prefix',metavar='str',action='store',type=str,nargs='?', help='specify the xmf file prefix')
	parser.add_argument('--repeat','-r',dest='repeat',action='store_const',const=True, help='use the first wedge for all slices')
	parser.add_argument('--quiet','-q',dest='quiet',action='store_const',const=True, help='only display error messages (default full debug messages)')
	parser.add_argument('--short','-s',dest='shortfilename',action='store_const',const=True, help='use shorter filenaming convention')
	parser.add_argument('--norepeat',dest='norepeat',action='store_const',const=True, help='debug variable for infinite recursive execution escaping')
	parser.add_argument('--disable',dest='disable',action='store',metavar='str',type=str, nargs='+', help='disable output of \'abundances\', \'hydro\',or \'radiation\' components AND/OR \'luminosity\' and/or \'E_RMS\' components')
	parser.add_argument('--xdmf',dest='xdmf',action='store_const',const=True, help='use .xdmf extension instead of default .xmf')
	parser.add_argument('--directory',dest='dir',metavar='str', const='.',action='store',type=str,nargs='?',help='Output xdmf in dirctory specified instead of next to hdf files')
	parser.add_argument('--auxiliary','-a',dest='aux',action='store_const',const=True, help='Write auxiliary computed (derivative) values like luminosity to a companion file')
	parser.add_argument('--reduce',dest='reduce',action='store_const',const=True, help='Also copy over the reduced data to the auxilary file')
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
			qprint("h5py load failed. Trying to run under reloaded modules")
			try:
				if not args.norepeat:
					if socket.gethostname()[:4]=='rhea':
						sp.call(['bash -cl "cd '+os.getcwd()+'; module unload PE-intel python;module load PE-gnu python python_h5py;python '+(' '.join(sys.argv))+' --norepeat"'],shell=True)
					if socket.gethostname()[:5]=='titan':
						sp.call(['bash -cl "cd '+os.getcwd()+'; module load python python_h5py;python '+(' '.join(sys.argv))+' --norepeat"'],shell=True)
					else:
						raise ValueError('aw crap, no known \'module\' command')
				else:
					raise ValueError('aw crap, already reloaded and still failed to import something')
			except:
				#check for module command
				sp.call(["module"],shell=True)
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
	try:
		mp.Semaphore()
		from joblib import Parallel, delayed
		num_cores=min(16,mp.cpu_count())
		if args.threads:
			num_cores=args.threads[0]
			qprint("Set threads to "+str(num_cores))
	except:
		eprint("Warning: joblib did not load correctly")
		qprint("	Attempting to override thread count to 1 and skip parallelization")
		num_cores=1
		if args.threads and args.threads!=1:
			eprint("	Warning: cannot override thread count, thread flag ignored")
	if args.repeat and num_cores!=1:
		num_cores=1
		qprint("Running with single thread")
	############################################################################################################################################################################################
	# On with bulk of code
	old_time=start_time # for speed diagnostics
	unique_set=set()
	for filename in args.files:
		unique_set.add(filename[:-5]+'01.h5')
	if len(unique_set)!=len(args.files):
		qprint('Warning: redundant files databases found and ignored \n\tOne or more foo_01.h5 and foo_0\d.h5 found in filename list')
	for filename in unique_set:
		hf = h5py.File(filename,'r')
		dims=[]
		for dim in hf['mesh']['array_dimensions']: #full dimensions of mesh
			dims.append(dim) 
		extents = []  #IE the dimensions actually used in this particular run
		extents.append(hf['mesh']['radial_index_bound'][1])
		extents.append(hf['mesh']['theta_index_bound'][1])
		extents.append(hf['mesh']['phi_index_bound'][1])
		topo_type='3DRectMesh'
		geom_type='VXVYVZ'
		is_3d=True
		is_2d=False
		if extents[2]==1:
			del extents[2]
			topo_type='2DRectMesh'
			geom_type=geom_type[0:4]
			is_2d=True
			is_3d=(not is_2d)
		n_hyperslabs = hf['mesh']['nz_hyperslabs'].value
		n_elemental_species = hf['abundance']['xn_c'].shape[3]
		entities={
			"Extent_r":str(extents[0]),
			"Extent_theta":str(extents[1]),
			"Dim_r":str(extents[0]+1),
			"h5path":str(filename[:-5])
		}
		if is_3d:
			entities['Extent_phi']=extents[2]
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
		extension='.xmf'
		if not file_directory.endswith('/') and file_directory != '':
			file_directory+='/'
		if args.xdmf:
			extension='.xdmf'
		if args.shortfilename:
			file_out_name=file_directory+filename_part[0]+'-'+filename_part[3]+'-'+filename_part[1]+extension
		else:
			file_out_name=file_directory+filename_part[0]+'_grid-'+filename_part[3]+'_step-'+filename_part[1]+extension
		# if input file doesn't have _pro suffix, rename input files to have this suffix. Later will write data if it is missing.
		processed_suffix='_pro'
		if filename[-7:-3]!='_pro':
			processed_suffix='' #temporary fix
		# 	old_filename=filename
		# 	filename=filename[0:-3]+'_pro.h5'
		# 	os.rename(old_filename,filename)
			# if not args.quiet:
				# print(old_filename+" renamed to "+filename+".")
				# print(Processed data added to h5 file)
		slices=n_hyperslabs #default slices value if not overidden by "--slices/-s int" argument on program call
		if args.slices:
			if not args.slices>slices:
				slices=args.slices
			else:
				eprint("Error: slices must not be more than the number of wedges")
				sys.exit()
		############################################################################################################################################################################################
		# compute luminosity auxilary variabes
		if args.aux or args.reduce:
			qprint("Creating auxillary file")
			aux_hf=h5py.File(re.sub("\d\d\.h5",'aux.h5',re.sub("\d\d_pro\.h5",'aux_pro.h5',filename)),'w')
			# prune info
			if args.reduce:
				def copygroup(group):
					br()
					sliced_quantities=[
					'abundance/*',
					'',
					'',
					]
					aux_hf.create_group(group[0])
					if slices>1:
						for item in hf[group[0]]:
							aux_hf.copy(group[0]+'/'+item,hf[group[0]][item])
							# if item in sliced_values:
							# 	for i in range(1,slices):
							# 		temp_hf = h5py.File(re.sub("\d\d\.h5",str(format(i, '02d'))+'.h5',re.sub("\d\d_pro\.h5",str(format(i, '02d'))+'_pro.h5',filename)),'r')


				for group in hf.items():
					copygroup(group)
			
		if args.aux:
			qprint("Creating derived values")
			if not args.disable or "radiation" not in args.disable:
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
				step=dims[1]/n_hyperslabs
				#open new auxilary hdf file or overite existing one. 
				
				aux_hf.create_group("/radiation") #or do nothing if exists
				######## Compute E_RMS_array (size N_species) of arrays (size N_groups) ##############
				# # initialize variables for parallel loop
				if not args.disable or "E_RMS" not in args.disable:
					psi0_c=hf['radiation']['psi0_c'] 
					def compute_E_RMS_array(sl):	
						sl+=1
						i=sl
						if args.repeat:
							i=1
						qprint("Computing E_RMS_[1.."+str(n_species)+"] for slice "+str(sl)+" from "+re.sub("\d\d\.h5",str(format(i, '02d'))+'.h5',re.sub("\d\d_pro\.h5",str(format(i, '02d'))+'_pro.h5',filename)))
						temp_hf= h5py.File(re.sub("\d\d\.h5",str(format(i, '02d'))+'.h5',re.sub("\d\d_pro\.h5",str(format(i, '02d'))+'_pro.h5',filename)),'r')
						psi0_c=temp_hf['radiation']['psi0_c'][:]
						row=np.empty((n_species,dims[2]/n_hyperslabs,dims[1],dims[0]))
						for n in range(0,n_species):
							numerator=np.sum(psi0_c[:,:,:,n]*e5de,axis=3)
							denominator=np.sum(psi0_c[:,:,:,n]*e3de,axis=3)
							row[n][:][:][:]=np.sqrt(numerator/(denominator+1e-100))
						return row
					if num_cores!=1:
						results = Parallel(n_jobs=num_cores)(delayed(compute_E_RMS_array)(sl) for sl in range(0,n_hyperslabs))
						#concatenate together the member of each array within E_RMS_ARRAY and write to auxilary HDF file
						qprint("Concatenating E_RMS_[0.."+str(n_species)+"] results...")
						E_RMS_array=np.hstack(results)
					else:
						E_RMS_array=np.empty((n_species,dims[2],dims[1],dims[0]))
						for sl in range(0,n_hyperslabs):
							E_RMS_array[:,sl*step:(sl+1)*step,:,:]=compute_E_RMS_array(sl)
					for n in range(0,n_species):
						qprint("Writing E_RMS_"+str(n)+" out to file")
						aux_hf.create_dataset("/radiation/E_RMS_"+str(n),data=E_RMS_array[n])
					del E_RMS_array
					if 'results' in locals():
						del results
				# ######## luminosity part ###########
				if not args.disable or "luminosity" not in args.disable:
					qprint("Computing luminosities")
					psi1_e=hf['radiation']['psi1_e']
					radius=hf['mesh']['x_ef'].value
					agr_e=hf['fluid']['agr_e'].value
					cell_area_GRcorrected=4*pi*radius**2/agr_e**4
					lumin=np.empty((dims[2],dims[1],dims[0]))
					lumin_array=[]
					def compute_luminosity(species,sl):
						j=sl
						if args.repeat:
							sl=0
						temp_hf= h5py.File(re.sub("\d\d\.h5",str(format(sl+1, '02d'))+'.h5',re.sub("\d\d_pro\.h5",str(format(sl+1, '02d'))+'_pro.h5',filename)),'r')
						psi1_e=temp_hf['radiation']['psi1_e']
						qprint("	On slice "+str(j+1)+" of "+str(n_hyperslabs)+" from "+re.sub("\d\d\.h5",str(format(sl+1, '02d'))+'.h5',re.sub("\d\d_pro\.h5",str(format(sl+1, '02d'))+'_pro.h5',filename)))
						del temp_hf
						return np.sum(psi1_e[:,:,:,species]*e3de, axis=3)*np.tile(cell_area_GRcorrected[1:dims[0]+1],(dims[2]/n_hyperslabs,dims[1],1))*(cvel*ecoef*1e-51)
					lumin_array=np.empty((n_species,dims[2],dims[1],dims[0]))
					for species in range(0,n_species):
						qprint("Computing luminosity for species "+str(species)+":")
						if num_cores!=1:
							lumin_array[species] = np.vstack(Parallel(n_jobs=num_cores)(delayed(compute_luminosity)(species,sl) for sl in range(species,n_hyperslabs)))
						else:
							for sl in range(0,n_hyperslabs):
								lumin_array[species,sl*step:(sl+1)*step,:,:]=compute_luminosity(0,sl)
					qprint("########################################")
					for n,lumin in enumerate(lumin_array):
						qprint("Writing luminosity species "+str(n)+" to auxilary hdf file")
						aux_hf.create_dataset("/radiation/Luminosity_"+str(n),data=lumin)
			# now stack mask for YY grid:
			qprint("########################################")
			qprint("Creating On_grid_mask")
			aux_hf.create_group("/mesh") #or do nothing if exists
			mask_slice = hf['mesh']['ongrid_mask'].value
			mask = np.dstack( mask_slice for i in range(0,extents[0]))
			del mask_slice
			qprint("Writing On_grid_mask to auxillary file")
			aux_hf.create_dataset("/mesh/mask",data=mask)
			aux_hf.close()
			del aux_hf
		############################################################################################################################################################################################
		dimstr_sub = str(dims[1])+" "+str(dims[0])
		extents_sub = str(extents[1])+" "+str(extents[0])
		if not is_2d:
			extents[2]=extents[2]/n_hyperslabs*slices
			dimstr_sub = str(dims[2]/n_hyperslabs)+" "+str(dims[1])+" "+str(dims[0])
			extents_sub = str(dims[2]/n_hyperslabs)+" "+str(extents[1])+" "+str(extents[0])
		dimstr = ' '.join([str(x) for x in dims[::-1]])
		extents_str = ' '.join([str(x) for x in extents[::-1]])
		# create xdmf element
		xdmf = et.Element("Xdmf", Version="2.0")
		# create Domain element
		domain = et.SubElement(xdmf,"Domain")

		############################################################################################################################################################################################
		#create main "Hydro" grid that will contain most scalars and Abundance grid that will contain n and p counts
		grid = {'Hydro':et.SubElement(domain,"Grid",Name="Hydro"),
				'Abundance':et.SubElement(domain,"Grid",Name="Abundance"),
				'Radiation':et.SubElement(domain,"Grid",Name="Radiation"),
				'Mesh':et.SubElement(domain,"Grid",Name="Mesh")}
		et.SubElement(grid['Hydro'],"Topology",TopologyType=topo_type,NumberOfElements=' '.join([str(x+1) for x in extents[::-1]]))
		geometry = et.SubElement(grid['Hydro'],"Geometry",GeometryType=geom_type)
		coords=["x_ef","y_ef","z_ef"]
		if is_2d:
			del coords[2]
		for n,coord_name in enumerate(coords):
			parent_element=geometry
			if coord_name=='x_ef':
				unit_changing_function = et.SubElement(geometry,"DataItem",Dimensions=str(extents[n]+1),ItemType="Function",Function="$0/100000")
				parent_element=unit_changing_function
			hyperslab = et.SubElement(parent_element, "DataItem",Dimensions=str(extents[n]+1),ItemType="HyperSlab")
			et.SubElement(hyperslab,"DataItem",Dimensions="3 1",Format="XML").text = "0 1 "+str(extents[n]+1)
			et.SubElement(hyperslab,"DataItem",Dimensions=str(hf['mesh'][coord_name].size),NumberType="Float",Precision="8",Format="HDF").text = "&h5path;01" + processed_suffix + ".h5:/mesh/" + coord_name
			
		et.SubElement(grid['Hydro'],"Time",Value=str(hf['mesh']['time'].value-hf['mesh']['t_bounce'].value))
		for name in grid:
			if name!='Hydro':
				et.SubElement(grid[name],"Topology",Reference="/Xdmf/Domain/Grid[1]/Topology[1]")
				et.SubElement(grid[name],"Geometry",Reference="/Xdmf/Domain/Grid[1]/Geometry[1]")

		############################################################################################################################################################################################
		# The following functions are helpers to create dimensions strings and "JOIN($0; $1; $2 .... $N<=9)" strings for the various hyperslabs and nested join functions
		m=dims[:]
		if not is_2d:
			m[2]=extents[2]
		block_string=' '.join([str(x) for x in m[::-1]])
		#block_string is a string that uses extents for the radial component and dims for the others
		def function_str(m):
			function_stri="JOIN("
			function_str_array=[]
			for n in range(0, int(m)):
				function_str_array.append("$"+str(n))
			function_stri+='; '.join(function_str_array)+")"
			return function_stri
		def dim_str(m):
			m=min([(slices-(m)*10),10])
			return str(dims[2]/n_hyperslabs*m)+" "+str(dims[1])+" "+str(dims[0])
		def extents_stri(m):
			return str(dims[2]/n_hyperslabs*m)+" "+str(extents[1])+" "+str(extents[0])
		############################################################################################################################################################################################
		# Loop through all standard scalars in "Hydro" grid
		at = et.SubElement(grid['Mesh'],"Attribute",Name="On_Grid_Mask",AttributeType="Scalar",Center="Cell",Dimensions=extents_str)
		et.SubElement(at,"DataItem",Dimensions=extents_str,NumberType="Float",Precision="8",Format="HDF").text= "&h5path;aux.h5:/mesh/mask"
		if not args.disable or "hydro" not in args.disable:
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
				"Density":"rho_c",
			}
			for name in storage_names:
				at = et.SubElement(grid['Hydro'],"Attribute",Name=name,AttributeType="Scalar",Center="Cell",Dimensions=extents_str)
				hyperslab = et.SubElement(at,"DataItem",Dimensions=extents_str,ItemType="HyperSlab")
				et.SubElement(hyperslab,"DataItem",Dimensions="3 3",Format="XML").text="0 0 0 1 1 1 "+[extents_str,'1 '+extents_str][is_2d]
				if is_3d:
					superfun = et.SubElement(hyperslab,"DataItem",ItemType="Function", Function=function_str(int((slices-1)/10)+1),Dimensions=block_string)
				n=1
				for m in range(0, int((slices-1)/10)+1):
					if is_3d:
						fun = et.SubElement(superfun,"DataItem",ItemType="Function", Function=function_str(min([(slices-m*10),10])),Dimensions=dim_str(m))
					for i in range(0,min(slices-m*10,10)):
						return_element=hyperslab
						if 'fun' in globals():
							return_element=fun
						if args.repeat:
							et.SubElement(return_element,"DataItem",Dimensions=[dimstr_sub,'1 '+dimstr_sub][is_2d],NumberType="Float",Precision="8",Format="HDF").text= "&h5path;01" + processed_suffix + ".h5:/fluid/" + storage_names[name]
						else:
							et.SubElement(return_element,"DataItem",Dimensions=[dimstr_sub,'1 '+dimstr_sub][is_2d],NumberType="Float",Precision="8",Format="HDF").text= "&h5path;" + str(format(n, '02d')) + processed_suffix + ".h5:/fluid/" + storage_names[name]
						n+=1
		############################################################################################################################################################################################
		# Abundance part:
		at = et.SubElement(grid['Abundance'],"Attribute",Name='nse_flag',AttributeType="Scalar",Center="Cell",Dimensions=extents_str)
		hyperslab = et.SubElement(at,"DataItem",Dimensions=extents_str,ItemType="HyperSlab")
		et.SubElement(hyperslab,"DataItem",Dimensions="3 3",Format="XML").text="0 0 0 1 1 1 "+[extents_str,'1 '+extents_str][is_2d]
		dim_nse=hf['abundance']['nse_c'].shape
		dim_nse_str=str(dim_nse[0])+' '+str(dim_nse[1])+' '+str(dim_nse[2])
		if is_3d:
			superfun = et.SubElement(hyperslab,"DataItem",ItemType="Function", Function=function_str(int((slices-1)/10)+1),Dimensions=str(dim_nse[0]*slices)+' '+str(dim_nse[1])+' '+str(dim_nse[2]))
		n=1
		for m in range(0, int((slices-1)/10)+1):
			if is_3d:
				fun = et.SubElement(superfun,"DataItem",ItemType="Function", Function=function_str(min([(slices-m),10])),Dimensions=str(dim_nse[0]*[slices%10,10][slices-m*10>=10])+" "+str(dim_nse[1])+" "+str(dim_nse[2]))
			for i in range(0,min(slices-m*10,10)):
				return_element=hyperslab
				if 'fun' in globals():
					return_element=fun
				if args.repeat:
					et.SubElement(return_element,"DataItem",Dimensions=dim_nse_str,NumberType="Int",Format="HDF").text= "&h5path;01" + processed_suffix + ".h5:/abundance/nse_c"
				else:
					et.SubElement(return_element,"DataItem",Dimensions=dim_nse_str,NumberType="Int",Format="HDF").text= "&h5path;" + str(format(n, '02d')) + processed_suffix + ".h5:/abundance/nse_c"
				n+=1
		del dim_nse,dim_nse_str

		species_names=hf['abundance']['a_name'].value
		if n_elemental_species-1==species_names.shape[0]:
			species_names=np.append(species_names,'aux')
		if not args.disable or "abundance" not in args.disable:
			for el,name in enumerate(species_names):
				if re.findall('\D\d',name): #if there is a transition between a non digit to a digit in the element name (IE in "li3" it would match because of the "i3")
					element_name=re.sub('\d','',name).capitalize() #set element_name to the capitalized element without the number
					name=re.sub('\D','',name) #find the transition between elements name and number
					if not 'Abundance'+'/'+element_name in grid: #If the grid for that element doesn't already exist, create it 
						grid['Abundance'+'/'+element_name]=et.SubElement(domain,"Grid",Name='Abundance'+'/'+element_name,GridType="Uniform")
						et.SubElement(grid['Abundance'+'/'+element_name],"Topology",Reference="/Xdmf/Domain/Grid[1]/Topology[1]")
						et.SubElement(grid['Abundance'+'/'+element_name],"Geometry",Reference="/Xdmf/Domain/Grid[1]/Geometry[1]")
					attribute=et.SubElement(grid['Abundance'+'/'+element_name],"Attribute",Name=name,AttributeType="Scalar",Center="Cell")
				else:
					attribute=et.SubElement(grid['Abundance'],"Attribute",Name=name,AttributeType="Scalar",Center="Cell")
				if is_3d:
					superfun = et.SubElement(attribute,"DataItem",ItemType="Function", Function=function_str(int((slices-1)/10)+1),Dimensions=extents_str)
				n=1
				for m in range(0, int((slices-1)/10)+1):
					if is_3d:
						fun = et.SubElement(superfun,"DataItem",ItemType="Function", Function=function_str(min([(slices-m*10),10])),Dimensions=extents_stri(min(slices-m*10,10)))
					for i in range(0,min(slices-m*10,10)):
						return_element=attribute
						if 'fun' in globals():
							return_element=fun
						dataElement = et.SubElement(return_element,"DataItem", ItemType="HyperSlab", Dimensions=extents_sub)
						et.SubElement(dataElement,"DataItem",Dimensions="3 4",Format="XML").text="0 0 0 "+str(el)+" 1 1 1 1 "+[extents_sub,'1 '+extents_str][is_2d]+" 1"
						if args.repeat==True:
							et.SubElement(dataElement,"DataItem",Dimensions=[dimstr_sub,dimstr][is_2d]+" "+str(n_elemental_species),NumberType="Float",Precision="8",Format="HDF").text= "&h5path;01" + processed_suffix + ".h5:/abundance/xn_c"
						else:
							et.SubElement(dataElement,"DataItem",Dimensions=[dimstr_sub,dimstr][is_2d]+" "+str(n_elemental_species),NumberType="Float",Precision="8",Format="HDF").text= "&h5path;" + str(format(n, '02d')) + processed_suffix + ".h5:/abundance/xn_c"
						n+=1
		############################################################################################################################################################################################
		##Create luminosity and E_RMS xdmf
		if not args.disable or "radiation" not in args.disable:
			n_species=hf['radiation']['raddim'][1] # in case the auxilary data is not generated this run
			for sp in ['e','e-bar','mt','mt-bar']:
				# E_RMS_[sp]:
				attribute = et.SubElement(grid['Radiation'],"Attribute",Name="E_RMS_"+sp,AttributeType="Scalar", Dimensions=extents_str,Center="Cell")
				hyperslab = et.SubElement(attribute, "DataItem",ItemType="HyperSlab",Dimensions=extents_str)
				et.SubElement(hyperslab,"DataItem",Dimensions="3 3",Format="XML").text="0 0 0 1 1 1 "+[extents_str,'1 '+extents_str][is_2d]
				et.SubElement(hyperslab,"DataItem",Dimensions=dimstr, Format="HDF").text= "&h5path;aux.h5:/radiation/E_RMS_"+sp
				# Luminosity_[sp]:
				attribute = et.SubElement(grid['Radiation'],"Attribute",Name="Luminosity_"+sp,AttributeType="Scalar", Dimensions=extents_str,Center="Cell")
				hyperslab = et.SubElement(attribute, "DataItem",ItemType="HyperSlab",Dimensions=extents_str)
				et.SubElement(hyperslab,"DataItem",Dimensions="3 3",Format="XML").text="0 0 0 1 1 1 "+[extents_str,'1 '+extents_str][is_2d]
				et.SubElement(hyperslab,"DataItem",Dimensions=dimstr, Format="HDF").text= "&h5path;aux.h5:/radiation/Luminosity_"+sp

		############################################################################################################################################################################################
		# Write document tree to file
		try:
			f=open(file_out_name,'w')
			del extension,file_directory
			# if lxml module loaded use it to write document (fasted, simplest implementation):
			entity_str = ''
			for key,value in six.iteritems(entities):
				entity_str+="\n  <!ENTITY "+key+" \""+str(value)+"\">"
			entity_str+="\n  <!-- Note that Dim_r must be exactly 1 more than Extent_r or VisIt will have a spontanious freak out session -->"
			try:
				#write to file:
				f.write(\
					#remove all the '&amp;' tags that appear due to xml parser and replace with '&' so the aliasing works
					str(\
						re.sub(\
							b'&amp;',b'&',et.tostring(xdmf,\
								pretty_print=True,\
								xml_declaration=True,\
								encoding="ASCII",\
								doctype="<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" ["+entity_str+"\n]>"\
								)\
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
				f.write("<?xml version='1.0' encoding='ASCII'?>\n<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" ["+entity_str+"\n]>")
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
		############################################################################################################################################################################################
	#end loop over all hdf5 files
	############################################################################################################################################################################################