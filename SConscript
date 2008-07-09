# -*- python -*-
# $Header$
# Authors: Leon Rochester <lsrea@slac.stanford.edu>
# Version: TkrUtil-03-14-00
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('TkrUtilLib', depsOnly = 1)
TkrUtil = libEnv.SharedLibrary('TkrUtil', listFiles(['src/*.cxx', 'src/Dll/*.cxx']))

progEnv.Tool('TkrUtilLib')
progEnv.Tool('EventLib')
test_TkrUtil = progEnv.GaudiProgram('test_TkrUtil', ['src/test/test_TkrUtil.cxx', 'src/test/test_TkrUtil_load.cxx'], test = 1)
test_IndexedVector = progEnv.Program('test_IndexedVector', 'src/test/testIndexedVector.cxx')

progEnv.Tool('registerObjects', package = 'TkrUtil', libraries = [TkrUtil], testApps = [test_TkrUtil, test_IndexedVector], includes = listFiles(['TkrUtil/*.h']))
