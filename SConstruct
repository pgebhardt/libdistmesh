import sconshelper

# create environment
env = Environment()

# use clang++
env.Replace(CXX='clang++')

# create librarty
sconshelper.Library(name='distmesh', env=env, arguments=ARGUMENTS,
    CXXFLAGS=[
        '-std=c++11',
        '-O4',
        ],
    LINKFLAGS=[
        '-O4',
        ],
    CPPPATH=[
        '/usr/include/eigen3/',
        ],
    LIBS=[
        'qhull',
        ],
    )
