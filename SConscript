# -*- python -*-
# $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/SConscript,v 1.25 2010/02/26 01:54:00 lsrea Exp $
# Authors: Leon Rochester <lsrea@slac.stanford.edu>
# Version: TkrUtil-03-19-02
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('TkrUtilLib', depsOnly = 1)
TkrUtil = libEnv.SharedLibrary('TkrUtil',
                               listFiles(['src/*.cxx', 'src/Dll/*.cxx']))

progEnv.Tool('TkrUtilLib')
progEnv.Tool('EventLib')
test_TkrUtil = progEnv.GaudiProgram('test_TkrUtil',
                                    ['src/test/test_TkrUtil.cxx',
                                     'src/test/test_TkrUtil_load.cxx'],
                                    test = 1)
test_IndexedVector = progEnv.Program('test_IndexedVector',
                                     ['src/test/testIndexedVector.cxx'])

progEnv.Tool('registerTargets', package = 'TkrUtil',
             libraryCxts = [[TkrUtil, libEnv]],
             testAppCxts = [[test_TkrUtil, progEnv],
                            [test_IndexedVector, progEnv]],
             includes = listFiles(['TkrUtil/*.h']))




