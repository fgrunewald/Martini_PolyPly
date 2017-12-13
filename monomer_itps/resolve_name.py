def resolve_name(names, path):
    paths = []
    for name in names:
        full_name = path + 'monomer_itps/' + name + '.itp'
        paths += [ full_name ]
    return(paths)

