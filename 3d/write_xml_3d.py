#!python
from __future__ import print_function
# For diagnostics, time the execution of the code
import time
start_time = time.time()
# Import all the things robustly
##############################################################################################
import sys, os, math, socket, argparse
import subprocess as sp
from pdb import set_trace as br

#define an error printing function for more accurate error reporting to terminal
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
# For ORNL
if socket.gethostname()[:4]=='rhea':
	sys.path.append('/lustre/atlas/proj-shared/ast109/amos/lib/python2.7/site-packages')
	sys.path.append('/sw/redhat6/visit/current/linux-x86_64/lib/site-packages/')
# For my laptop
elif socket.gethostname()=='Lycoris':
	sys.path.append('/Applications/VisIt.app/Contents/Resources/2.10.2/darwin-x86_64/lib/site-packages')
try:
	import visit
	import visit_utils
except ImportError:
	eprint("Error: visit module import failed\n")
	eprint('Please make sure visit\'s path is in $PATH')
	eprint('	( /path/to/visit/bin )')
	eprint('Please make sure visit\'s /lib/site-packages directory is in $PYTHONPATH')
	eprint ('	( /path/to/visit/VersionNumber/platform/lib/site-packages )')
	sys.exit()
# Contstruct Parser
##############################################################################################
parser = argparse.ArgumentParser(description="Generate XDMF files from Chimera h5 files")
# parser.files is a list of 1 or more h5 files that will have xdmf files generated for them
parser.add_argument('files',metavar='foo.h5',type=str,nargs='+',help='h5 files to process (1 or more args)')
parser.add_argument('--extents','-e',dest='dimensions',metavar='int',action='store',type=int,nargs=3, help='dimensions to crop to')
parser.add_argument('--slices','-s',dest='slices',metavar='int',action='store',type=int,nargs=1, help='number of slices to use')
parser.add_argument('--prefix','-p',dest='prefix',metavar='str',action='store',type=str,nargs=1, help='specify the xmf file prefix')
parser.add_argument('--repeat','-r',dest='repeat',action='store_const',const=True, help='use the first wedge for all slices')
parser.add_argument('--quiet','-q',dest='quiet',action='store_const',const=True, help='use the first wedge for all slices')
parser.add_argument('--disable',dest='disable',action='store_const',const=True, help='debug variable to make script not run again on failure')
args=parser.parse_args()
##############################################################################################
#This next bit is specific to ORNL. If h5py import fails it switches environments and reloads this script
try:
	#Most likly to fail part.
	import h5py
# Try and correct the h5import error by launching subprocess that calls this script again after loading proper modules on rhea
except ImportError:
	try:
		if not args.quiet:
			print("Trying to run under reloaded modules")
		try:
			if not args.run:
				sp.call(["module unload PE-intel python;module load PE-gnu python python_h5py"],shell=True)
			else:
				raise ValueError('aw crap')
		except:
			#redo the offending call so the error can display
			sp.call(["module unload PE-intel python;module load PE-gnu python python_h5py"],shell=True)
			eprint("Could not import modules")
			raise ValueError('aw crap')
		sp.call(["module unload PE-intel python;module load PE-gnu python python_h5py;python "+(' '.join(sys.argv))+" --run"],shell=True)
		if not args.quiet:
			print("Finished")
	except:
		eprint("Fatal error: could not import h5py or reload modules to make it possible. h5 reading is impossible without h5py.")
	sys.exit()
#Robustly import an xml writer/parser
try:
	from lxml import etree as et
	if not args.quiet:
		print("Running with lxml.etree")
except ImportError:
	try:
		# Python 2.5
		import xml.etree.cElementTree as et
		import xml.dom.minidom as md
		if not args.quiet:
			print("Running with cElementTree on Python 2.5+")
	except ImportError:
		try:
		# Python 2.5
			import xml.etree.ElementTree as et
			if not args.quiet:
				print("Running with ElementTree on Python 2.5+")
		except ImportError:
			try:
				# normal cElementTree install
				import cElementTree as et
				if not args.quiet:
					print("running with cElementTree")
			except ImportError:
				try:
					# normal ElementTree install
					import elementtree.ElementTree as et
					if not args.quiet:
						print("running with ElementTree")
				except ImportError:
					eprint("Fatal error: Failed to import ElementTree from any known place. XML writing is impossible. ")

import numpy as np
import re
#End Parser construction
# On with bulk of code
old_time=start_time
for filename in args.files:
	hf = h5py.File(filename,'r')
	dims=[]
	for dim in hf['mesh']['array_dimensions']: #full dimensions of mesh
		dims.append(dim) 
	extents = []  #IE the dimensions actually used in this particular run
	extents.append(hf['mesh']['radial_index_bound'][1])
	extents.append(hf['mesh']['theta_index_bound'][1])
	extents.append(hf['mesh']['phi_index_bound'][1])
	n_hyperslabs = hf['mesh']['nz_hyperslabs'].value

	slices=n_hyperslabs #default slices value if not overidden by "--slices/-s int" argument on program call
	if args.slices:
		if not args.slices[0]>slices:
			slices=args.slices[0]
		else:
			eprint("Error: slices must not be more than the number of wedges")
			sys.exit()
	extents[2]=extents[2]/n_hyperslabs*slices
	dimstr_sub = str(dims[2]/n_hyperslabs)+" "+str(dims[1])+" "+str(dims[0])
	extents_sub = str(dims[2]/n_hyperslabs)+" "+str(extents[1])+" "+str(extents[0])
	dimstr = ' '.join([str(x) for x in dims[::-1]])
	extents_str = ' '.join([str(x) for x in extents[::-1]])
	# create xdmf element
	xdmf = et.Element("Xdmf", Version="2.0")
	# create Domain element
	domain = et.SubElement(xdmf,"Domain")

	#create main "Hydro" grid that will contain most scalars and Abundance grid that will contain n and p counts 
	grid = {'Hydro':et.SubElement(domain,"Grid",Name="Hydro"),
			'Abundance':et.SubElement(domain,"Grid",Name="Abundance")}
	et.SubElement(grid['Hydro'],"Topology",TopologyType="3DRectMesh",NumberOfElements=' '.join([str(x+1) for x in extents[::-1]]))
	geometry = et.SubElement(grid['Hydro'],"Geometry",GeometryType="VXVYVZ")
	coords=["x_ef","y_ef","z_ef"]
	for n,coord_name in enumerate(coords):
		hyperslab = et.SubElement(geometry, "DataItem",Dimensions=str(extents[n]+1),ItemType="HyperSlab",Type="HyperSlab")
		et.SubElement(hyperslab,"DataItem",Dimensions="3 1",Format="XML").text="0 1 "+str(extents[n]+1)
		et.SubElement(hyperslab,"DataItem",Dimensions=str(hf['mesh'][coord_name].size),NumberType="Float",Precision="8",Format="HDF").text = filename + ":/mesh/" + coord_name
		
	et.SubElement(grid['Hydro'],"Time",Value=str(hf['mesh']['time'].value-hf['mesh']['t_bounce'].value))
	et.SubElement(grid['Abundance'],"Topology",Reference="/Xdmf/Domain/Grid[1]/Topology[1]")
	et.SubElement(grid['Abundance'],"Geometry",Reference="/Xdmf/Domain/Grid[1]/Geometry[1]")
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
		# "mean_A":"",#computed quantity to be added later
		"nse_flag":"e_int",
	}
	# The following is functions are helpers to create dimensions strings and "JOIN($0; $1; $2 .... $N<=9)" strings for the various hyperslabs and nested join functions
	m=dims[:]
	m[2]=extents[2]
	block_string=' '.join([str(x) for x in m[::-1]])
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
	# Loop through all standard scalars in "Hydro" grid
	for name in storage_names:
		at = et.SubElement(grid['Hydro'],"Attribute",Name=name,AttributeType="Scalar",Center="Cell",Dimensions=extents_str)
		hyperslab = et.SubElement(at,"DataItem",Dimensions=extents_str,ItemType="HyperSlab",Type="HyperSlab")
		et.SubElement(hyperslab,"DataItem",Dimensions="3 3",Format="XML").text="0 0 0 1 1 1 "+extents_str
		superfun = et.SubElement(hyperslab,"DataItem",ItemType="Function", Function=function_str(int((slices-1)/10)+1),Dimensions=block_string)
		n=1
		for m in range(0, int((slices-1)/10)+1):
			fun = et.SubElement(superfun,"DataItem",ItemType="Function", Function=function_str(min([(slices-m*10),10])),Dimensions=dim_str(m))
			for i in range(0,min(slices-m*10,10)):
				if args.repeat:
					et.SubElement(fun,"DataItem",Dimensions=dimstr_sub,NumberType="Float",Precision="8",Format="HDF").text= filename + ":/fluid/" + storage_names[name]
				else:
					et.SubElement(fun,"DataItem",Dimensions=dimstr_sub,NumberType="Float",Precision="8",Format="HDF").text= filename[:-5] + str(format(n, '02d')) + ".h5:/fluid/" + storage_names[name]
				n+=1	
	# Now loop through all the abundance elements and generate hyperslabs
	for el,name in enumerate(hf['abundance']['a_name']):
		if re.findall('\D\d',name):
			element_name=re.sub('\d','',name)
			name=re.sub('\D','',name) #find the transition between elements name and number
			#If the grid for that element doesn't already exist create it 
			if not grid.has_key('Abundance'+'/'+element_name):
				grid['Abundance'+'/'+element_name]=et.SubElement(domain,"Grid",Name='Abundance'+'/'+element_name,GridType="Uniform")
				et.SubElement(grid['Abundance'+'/'+element_name],"Topology",Reference="/Xdmf/Domain/Grid[1]/Topology[1]")
				et.SubElement(grid['Abundance'+'/'+element_name],"Geometry",Reference="/Xdmf/Domain/Grid[1]/Geometry[1]")
			attribute=et.SubElement(grid['Abundance'+'/'+element_name],"Attribute",Name=name,AttributeType="Scalar",Center="Cell")
		else:
			attribute=et.SubElement(grid['Abundance'],"Attribute",Name=name,AttributeType="Scalar",Center="Cell")
		superfun = et.SubElement(attribute,"DataItem",ItemType="Function", Function=function_str(int((slices-1)/10)+1),Dimensions=extents_str)
		n=1
		for m in range(0, int((slices-1)/10)+1):
			fun = et.SubElement(superfun,"DataItem",ItemType="Function", Function=function_str(min([(slices-m*10),10])),Dimensions=extents_stri(min(slices-m*10,10)))
			for i in range(0,(slices-m*10)):
				dataElement = et.SubElement(fun,"DataItem", ItemType="HyperSlab", Dimensions=extents_sub, Type="HyperSlab")
				et.SubElement(dataElement,"DataItem",Dimensions="3 4",Format="XML").text="0 0 0 "+str(el)+" 1 1 1 1 "+extents_sub+" 1"
				if args.repeat==True:
					et.SubElement(dataElement,"DataItem",Dimensions=dimstr_sub+" 17",NumberType="Float",Precision="8",Format="HDF").text= filename + ":/abundance/xn_c"
				else:
					et.SubElement(dataElement,"DataItem",Dimensions=dimstr_sub+" 17",NumberType="Float",Precision="8",Format="HDF").text= filename[:-5] + str(format(n, '02d')) + ".h5:/abundance/xn_c"
				n+=1

	# Write document tree to file
	try:
		f=open(filename[:-6]+'.xmf','w')
		# if lxml module loaded
		try:
			f.write(et.tostring(xdmf,pretty_print=True,xml_declaration=True,doctype="<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>"))
		#other ElementTree writers
		except:
			f.close()
			f=open(filename[:-6]+'.xmf','w')
			if not args.quiet:
				print("Writing "+filename+".xmf with improvised \"pretty print\"")
			def prettify(elem):
				rough_string = et.tostring(elem, 'utf-8')
				reparsed = md.parseString(rough_string)
				t = reparsed.toprettyxml(indent="\t")
				t = ''.join(t.splitlines(True)[1:]) #removes extra doc declaration
				return t
			# write custom doctype declaration
			f.write("<?xml version='1.0' encoding='ASCII'?>\n<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n")
			f.close()
			f=open(filename[:-6]+'.xmf','a')
			f.write(prettify(xdmf))
		f.close()
		if not args.quiet:
			print("--- XMF file created in %s seconds ---" % (time.time()-old_time))
		old_time=time.time()
	except:
		eprint("Fatal error:")
		sys.exit(["	File write error. Try adding adding a br() (local name of pdb.set_trace()) in the writer section at the end of script to debug"])