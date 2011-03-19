# -*- python -*-
# $Header: /nfs/slac/g/glast/ground/cvs/TkrUtil/SConscript,v 1.42 2010/12/14 20:37:29 lsrea Exp $
# Authors: Leon Rochester <lsrea@slac.stanford.edu>
# Version: TkrUtil-03-21-08
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('addLinkDeps', package='TkrUtil', toBuild='component')
TkrUtil = libEnv.SharedLibrary('TkrUtil',
                               listFiles(['src/*.cxx', 'src/Dll/*.cxx']))

progEnv.Tool('TkrUtilLib')
progEnv.Tool('EventLib')
test_TkrUtil = progEnv.GaudiProgram('test_TkrUtil',
                                    ['src/test/test_TkrUtil.cxx',
                                     'src/test/test_TkrUtil_load.cxx'],
                                    test = 1, package='TkrUtil')
test_IndexedVector = progEnv.Program('test_IndexedVector',
                                     ['src/test/testIndexedVector.cxx'])

progEnv.Tool('registerTargets', package = 'TkrUtil',
             libraryCxts = [[TkrUtil, libEnv]],
             testAppCxts = [[test_TkrUtil, progEnv],
                            [test_IndexedVector, progEnv]],
             includes = listFiles(['TkrUtil/*.h']),
             jo = listFiles(['src/test/*.txt', 'src/test/*.xml']))




