# $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrUtil/TkrUtilLib.py,v 1.1 2008/07/09 21:13:43 glastrm Exp $
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['TkrUtil'])
    env.Tool('CalibDataLib')
    env.Tool('xmlBaseLib')
    env.Tool('geometryLib')
    env.Tool('EventLib')

def exists(env):
    return 1;
