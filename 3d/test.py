import sys
paths=['/Applications/VisIt.app/Contents/Resources/2.10.2/darwin-x86_64/lib/python/lib/python2.7',
'/Applications/VisIt.app/Contents/Resources/2.10.2/darwin-x86_64/lib/python/lib/python2.7/lib-dynload',
'/Applications/VisIt.app/Contents/Resources/2.10.2/darwin-x86_64/lib/python/lib/python27.zip',
'/Applications/VisIt.app/Contents/Resources/2.10.2/darwin-x86_64/lib/python/lib/python2.7/plat-darwin',
'/Applications/VisIt.app/Contents/Resources/2.10.2/darwin-x86_64/lib/python/lib/python2.7/plat-mac',
'/Applications/VisIt.app/Contents/Resources/2.10.2/darwin-x86_64/lib/python/lib/python2.7/plat-mac/lib-scriptpackages',
'/Applications/VisIt.app/Contents/Resources/2.10.2/darwin-x86_64/lib/python/lib/python2.7/lib-tk',
'/Applications/VisIt.app/Contents/Resources/2.10.2/darwin-x86_64/lib/python/lib/python2.7/lib-old',
'/Applications/VisIt.app/Contents/Resources/2.10.2/darwin-x86_64/lib/python/lib/python2.7/site-packages',
'/Applications/VisIt.app/Contents/Resources/2.10.2/darwin-x86_64/lib/python/lib/python2.7/site-packages/PIL',
'/Applications/VisIt.app/Contents/Resources/2.10.2/darwin-x86_64/lib/site-packages']
for path in paths:
	sys.path.append(path)
