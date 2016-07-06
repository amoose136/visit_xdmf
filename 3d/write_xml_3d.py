import time
start_time = time.time()
import sys
import os
import subprocess as sp
visit_path=sp.check_output('which visit',shell=True).rstrip()[:-9]
sys.path.append('/sw/redhat6/visit/current/linux-x86_64/lib/site-packages/')
sys.path.append('/Applications/VisIt.app/Contents/Resources/2.10.2/darwin-x86_64/lib/site-packages')
os.environ['DYLD_LIBRARY_PATH']
try:
	import visit
	import visit_utils
except ImportError:
	print("Error: visit module import faile\n")
	print('Please make sure visit\'s path is in $PATH')
	print('	( /path/to/visit/bin )')
	print('Please make sure visit\'s /lib/site-packages directory is in $PYTHONPATH')
	print ('	( /path/to/visit/VersionNumber/platform/lib/site-packages )')
	sys.exit()
visit.Launch()
try:
	filename=sys.argv[1]
except:
	print("Error: No filename argument provided.")
	print("Try instead:")
	print("	python write_xml.py foo.h5")
	sys.exit()
#This next bit is specific to ORNL. If h5py import fails it switches environments and reloads this script
try:
	#Most likly to fail part. On failure, go to exception.
	import h5py
	#Robustly import an xml writer/parser
	try:
		from lxml import etree as et
		print("Running with lxml.etree")
	except ImportError:
		try:
			# Python 2.5
			import xml.etree.cElementTree as et
			import xml.dom.minidom as md
			print("Running with cElementTree on Python 2.5+")
		except ImportError:
			try:
			# Python 2.5
				import xml.etree.ElementTree as et
				print("Running with ElementTree on Python 2.5+")
			except ImportError:
				try:
					# normal cElementTree install
					import cElementTree as et
					print("running with cElementTree")
				except ImportError:
					try:
						# normal ElementTree install
						import elementtree.ElementTree as et
						print("running with ElementTree")
					except ImportError:
						print("Failed to import ElementTree from any known place")

	import numpy as np
	import re

	# On with bulk of code

	hf = h5py.File(filename,'r')
	dims=[]

	#appends the first element of 'array_dimensions' -> 542 (722 in 3d set)
	dims.append(hf['mesh']['array_dimensions'][0]) #dims[0]=722
	#appends the second element of 'array_dimensions' -> 180 (8 in 3d set)
	dims.append(hf['mesh']['array_dimensions'][1]) #dims[1]=8
	#appends the third element of 'array_dimensions' -> 1 (16 in 3d set)
	dims.append(hf['mesh']['array_dimensions'][2]) #dims[2]=16
	dimstr = str(dims[2])+" "+str(dims[1])+" "+str(dims[0])
	dimstr_sub = str(dims[2]/hf['mesh']['nz_hyperslabs'].value)+" "+str(dims[1])+" "+str(dims[0])#dimstr=""
	# create xdmf element

	# create Domain element
	xdmf = et.Element("Xdmf",Version="2.0")
	domain = et.SubElement(xdmf,"Domain")

	grid = {'Hydro':et.SubElement(domain,"Grid",Name="Hydro"),
			'Abundance':et.SubElement(domain,"Grid",Name="Abundance")}
	et.SubElement(grid['Hydro'],"Topology",TopologyType="3DRectMesh",NumberOfElements=str(dims[2]+1)+" "+str(dims[1]+1)+" "+str(dims[0]+1))
	geometry = et.SubElement(grid['Hydro'],"Geometry",GeometryType="VXVYVZ")
	et.SubElement(geometry,"DataItem",Dimensions=str(dims[0]+1),NumberType="Float",Precision="8",Format="HDF").text = filename + ":/mesh/x_ef"
	et.SubElement(geometry,"DataItem",Dimensions=str(dims[1]+1),NumberType="Float",Precision="8",Format="HDF").text = filename + ":/mesh/y_ef"
	et.SubElement(geometry,"DataItem",Dimensions=str(dims[2]+1),NumberType="Float",Precision="8",Format="HDF").text = filename + ":/mesh/z_ef"
	et.SubElement(grid['Hydro'],"Time",Value=str(hf['mesh']['time'].value-hf['mesh']['t_bounce'].value))
	et.SubElement(grid['Abundance'],"Topology",Reference="/Xdmf/Domain/Grid[1]/Topology[1]")
	et.SubElement(grid['Abundance'],"Geometry",Reference="/Xdmf/Domain/Grid[1]/Geometry[1]")

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
	n_hyperslabs = hf['mesh']['nz_hyperslabs'].value
	function_str="JOIN("
	for n in range(0, n_hyperslabs):
			function_str+="$"+str(n)
			if n!=n_hyperslabs-1:
				function_str+=" ; "
			else:
				function_str+=")"
	for name in storage_names:
		at = et.SubElement(grid['Hydro'],"Attribute",Name=name,AttributeType="Scalar",Center="Cell",Dimensions=dimstr)
		fun = et.SubElement(at,"DataItem",ItemType="Function", Function=function_str,Dimensions=dimstr)
		for n in range(1, n_hyperslabs+1):
			et.SubElement(fun,"DataItem",Dimensions=dimstr_sub,NumberType="Float",Precision="8",Format="HDF").text= filename[:-5] + str(format(n, '02d')) + ".h5:/fluid/" + storage_names[name]

	for i,name in enumerate(hf['abundance']['a_name']):
		if re.findall('\D\d',name):
			element_name=re.sub('\d','',name)
			name=re.sub('\D','',name) #find the transition between elements name and number
			if not grid.has_key('Abundance'+'/'+element_name):
				grid['Abundance'+'/'+element_name]=et.SubElement(domain,"Grid",Name='Abundance'+'/'+element_name,GridType="Uniform")
				et.SubElement(grid['Abundance'+'/'+element_name],"Topology",Reference="/Xdmf/Domain/Grid[1]/Topology[1]")
				et.SubElement(grid['Abundance'+'/'+element_name],"Geometry",Reference="/Xdmf/Domain/Grid[1]/Geometry[1]")
			attribute=et.SubElement(grid['Abundance'+'/'+element_name],"Attribute",Name=name,AttributeType="Scalar",Center="Cell")
		else:
			attribute=et.SubElement(grid['Abundance'],"Attribute",Name=name,AttributeType="Scalar",Center="Cell")
		fun = et.SubElement(attribute,"DataItem",ItemType="Function", Function=function_str,Dimensions=dimstr)
		for n in range(1, n_hyperslabs+1):
			dataElement = et.SubElement(fun,"DataItem", ItemType="HyperSlab", Dimensions=dimstr_sub, Type="HyperSlab")
			et.SubElement(dataElement,"DataItem",Dimensions="3 4",Format="XML").text="0 0 0 "+str(i)+" 1 1 1 1 "+dimstr_sub+" 1"
			et.SubElement(dataElement,"DataItem",Dimensions=dimstr_sub+" 17",Precisions="8",Format="HDF").text=filename[:-5] + str(format(n, '02d'))+".h5:/abundance/xn_c"
	
	# Write document tree to file
	f=open(filename[:-3]+'.xmf','w')
	# if lxml
	try:
		f.write(et.tostring(xdmf,pretty_print=True,xml_declaration=True,doctype="<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>"))
	#others
	except:
		print("Writing "+filename+".xmf with improvised \"pretty print\"")
		def prettify(elem):
			rough_string = et.tostring(elem, 'utf-8')
			reparsed = md.parseString(rough_string)
			t = reparsed.toprettyxml(indent="\t")
			return t
		f.write("<?xml version='1.0' encoding='ASCII'?>\n<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n")
		f.close()
		f=open(filename[:-3]+'.xmf','a')
		f.write(prettify(xdmf))
	f.close()
	print("--- XMF file created in %s seconds ---" % (time.time()-start_time))
except ImportError:
	try:
		print("Trying to run under reloaded modules")
		try:
			sp.call(["module unload PE-intel python;module load PE-gnu python python_h5py"],shell=True)
			sp.call(["module unload PE-intel python;module load PE-gnu python python_h5py;python write_xml.py "+filename],shell=True)
		except:
			#redo the offending call so the error can display
			sp.call(["module unload PE-intel python;module load PE-gnu python python_h5py"],shell=True)
			print("Could not import modules")
		print("Finished")
	except:
		print("Fatal error: could not import h5py")
