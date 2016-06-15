import time
start_time =time.time()
#Taken from the lxml tutorial in the lxml manual because it is a robust import statement
try:
	from lxml import etree as et
	print("running with lxml.etree")
except ImportError:
	try:
		# Python 2.5
		import xml.etree.cElementTree as et
		print("running with cElementTree on Python 2.5+")
	except ImportError:
		try:
			# Python 2.5
			import xml.etree.ElementTree as et
			print("running with ElementTree on Python 2.5+")
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
# End directly coppied portion

import h5py
import numpy as np



filename='chimera_00774_grid_1_01.h5'
hf = h5py.File(filename,'r')
dims=[]
#appends the first element of 'array_dimensions' -> 542
dims.append(str(hf['mesh']['array_dimensions'][0])) #dims[0]=542
#appends the first element of 'array_dimensions' -> 180
dims.append(str(hf['mesh']['array_dimensions'][1])) #dims[1]=180
dimstr = str(dims[1])+" "+str(dims[0]) #dimstr="180 542"
grid = et.Element("Grid",Name="Abundance",GridType="Uniform")
et.SubElement(grid,"Topology",Reference="/XDMF/Domain/Grid[1]/Topology[1]")
et.SubElement(grid,"Grid",Reference="/XDMF/Domain/Grid[1]/Geometry[1]")
for i,name in enumerate(hf['abundance']['a_name'][:5]):
	attribute=et.SubElement(grid,"Attribute",Name=name,AttributeType="Scalar",Center="Cell")
	dataElement = et.SubElement(attribute,"DataItem", ItemType="HyperSlab", Dimensions=dimstr, Type="HyperSlab")
	et.SubElement(dataElement,"DataItem",Dimensions="3 4",Format="XML").text="\n0 0 0 "+str(i)+"\n1 1 1 1\n1 "+dimstr+" 1\n"
	et.SubElement(dataElement,"DataItem",Dimensions="1 "+dimstr+" 161",Precisions="8",Format="HDF").text="\n"+filename+":/abundance/xn_c\n"
f=open(filename[:-3]+'.xmf','w')
f.write(et.tostring(grid,pretty_print=True,xml_declaration=True))
print("--- XMF file created in %s seconds ---" % (time.time()-start_time))
