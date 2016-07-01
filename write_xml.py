import time
start_time =time.time()
import sys
try:
	filename=sys.argv[1]
except:
	print("Error: No filename argument provided.")
	print("Try instead:")
	print("	python write_xml.py foo.h5")
	sys.exit()
#This next bit is specific to ORNL. If h5py import fails it switches environments and reloads the file
try:
	#Most likly to fail part
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
	#appends the first element of 'array_dimensions' -> 542
	dims.append(str(hf['mesh']['array_dimensions'][0])) #dims[0]=542
	#appends the first element of 'array_dimensions' -> 180
	dims.append(str(hf['mesh']['array_dimensions'][1])) #dims[1]=180
	dimstr = str(dims[1])+" "+str(dims[0]) #dimstr="180 542"
	# create xdmf element

	# create Domain element
	xdmf = et.Element("Xdmf",Version="2.0")
	domain = et.SubElement(xdmf,"Domain")

	grid = {'Hydro':et.SubElement(domain,"Grid",Name="Hydro"),
			'Abundance':et.SubElement(domain,"Grid",Name="Abundance")}
	et.SubElement(grid['Hydro'],"Topology",TopologyType="2DRectMesh",NumberOfElements="181 543")
	geometry = et.SubElement(grid['Hydro'],"Geometry",GeometryType="VXVY")
	et.SubElement(geometry,"DataItem",Dimensions="543",NumberType="Float",Precision="8",Format="HDF").text = filename + ":/mesh/x_ef"
	et.SubElement(geometry,"DataItem",Dimensions="181",NumberType="Float",Precision="8",Format="HDF").text = filename + ":/mesh/y_ef"
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
	for name in storage_names:
		at = et.SubElement(grid['Hydro'],"Attribute",Name=name,AttributeType="Scalar",Center="Cell")
		et.SubElement(at,"DataItem",Dimensions=dimstr,NumberType="Float",Precision="8",Format="HDF").text= filename + ":/fluid/" + storage_names[name]

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
		dataElement = et.SubElement(attribute,"DataItem", ItemType="HyperSlab", Dimensions=dimstr, Type="HyperSlab")
		et.SubElement(dataElement,"DataItem",Dimensions="3 4",Format="XML").text="0 0 0 "+str(i)+" 1 1 1 1 1 "+dimstr+" 1"
		et.SubElement(dataElement,"DataItem",Dimensions="1 "+dimstr+" 161",Precisions="8",Format="HDF").text=filename+":/abundance/xn_c"
	
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
		import subprocess as sp
		print("Trying to run under reloaded modules")
		try:
			sp.call(["module unload PE-intel python;module load PE-gnu python python_h5py"],shell=True)
			sp.call(["module unload PE-intel python;module load PE-gnu python python_h5py;python write_xml.py "+filename],shell=True)
		except:
			sp.call(["module unload PE-intel python;module load PE-gnu python python_h5py"],shell=True)
			print("Could not import modules")
		print("Finished")
	except:
		print("Fatal error: could not import h5py")
