from __future__ import print_function
# global sys, os, math, socket, argparse, re, sp, np, br, mp, h5py, six, et
try:
	import builtins
	args=builtins.args
except: 
	import __builtin__
	args=__builtin__.args

def eprint(*arg, **kwarg):
	print(*arg, file=sys.stderr, **kwarg)
# #define a standard printing function that only functions if there is no silence flag on script invocation
def qprint(*arg,**kwargs):
	if not args.quiet:
		print(*arg,**kwargs)

import sys, os, math, socket, argparse, re
import subprocess as sp
import numpy as np
from pdb import set_trace as br #For debugging I prefer the c style "break" nomenclature to "trace"
import multiprocessing as mp #For parallel speedup in derivative values 
# For ORNL
if socket.gethostname()[:4]=='rhea':
	# sys.path.append('/lustre/atlas/proj-shared/ast109/amos/lib/python2.7/site-packages')
	sys.path.append('./lib/python2.7/site-packages') #added for the assumption that script is called from bas directory
	sys.path.append('/sw/redhat6/visit/current/linux-x86_64/lib/site-packages/')
if socket.gethostname()[:5]=='titan':
	sys.path.append('/lustre/atlas/proj-shared/ast109/amos/lib/python2.7/site-packages')
# For my laptop
elif socket.gethostname()=='Lycoris':
	sys.path.append('/Applications/VisIt.app/Contents/Resources/2.10.2/darwin-x86_64/lib/site-packages')
import six
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