# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['TkrUtil'])
    env.Tool('CalibDataLib')
    env.Tool('xmlBaseLib')

def exists(env):
    return 1;
