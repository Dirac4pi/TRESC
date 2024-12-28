'''
script for visualizing 2-component complex MOs

input: 2 Gaussian cube format files contain 
(alpha real & beta real), (alpha imagine & 
beta imagine)of a selected 2-component com-
plex MO

output: amplitude isosurface with phase colour
mapping, contains alpha & beta components
'''

from os import environ
import numpy as np
import qcelemental as qcel
qcel.CovalentRadii("ALVAREZ2008")
from json import load as jsload
from mayavi import mlab


def read_cube(cube_address, imo=0):
    '''
    load data in cube format file.
    ----------
    cube_address : TYPE string
    ---address of .cube file.
    imo : TYPE integer
    ---MO serial number to be read.
    ----------
    Returns np.array*2, integer*1
    '''
    try:
        f = open(cube_address, 'r')
    except FileNotFoundError:
        if cube_address.endswith('.cub'):
            cube_address = cube_address + 'e'
        elif cube_address.endswith('.cube'):
            cube_address = cube_address.rstrip('e')
    else:
        f.close()
    with open(cube_address, 'r') as f:
        # read natom and coordination
        while True:
            try:
                ncenter, orgx, orgy, orgz = map(float, f.readline().split())
            except ValueError:
                continue
            else:
                break
                
        # read grid information
        n1, v1x, v1y, v1z = map(float, f.readline().split())
        n2, v2x, v2y, v2z = map(float, f.readline().split())
        n3, v3x, v3y, v3z = map(float, f.readline().split())
        
        # determine if it contains MOs
        nmo = 0
        if ncenter < 0:
            nmo = 1
            ncenter = abs(ncenter)  # reduce to positive
        
        # allocate memory
        atoms = np.zeros(int(ncenter), dtype=[('index','i4'), \
            ('charge','f8'), ('x','f8'), ('y','f8'), ('z','f8')])
        cubmat = np.zeros((int(n1),int(n2),int(n3)), dtype=[('value','f8'), \
            ('x','f8'), ('y','f8'), ('z','f8')])

        # read atomic information
        for i in range(int(ncenter)):
            atoms[i]['index'], atoms[i]['charge'], atoms[i]['x'], \
                atoms[i]['y'], atoms[i]['z'] = map(float, f.readline().split())

        if nmo == 1:  # read the number of MOs if exist
            nmo = int(f.readline().split()[0])
            if nmo < imo:
                exit('bad call read_cube: nmo < imo.')
            elif imo <= 0:
                exit('bad call read_cube: nmo > 0 but imo <= 0.')
            # read grid data
            for i in range(int(n1)):
                for j in range(int(n2)):
                    iline = 1
                    line = f.readline().split()
                    for k in range(imo, int(n3*nmo)+1, nmo):
                        if k // 6 >= iline:
                            line = f.readline().split()
                            iline += 1
                        if imo == nmo:
                            cubmat[i,j,k//nmo-1]['value'] = float(line[k%6-1])
                            #line[-1] is same as line[5]
                        else:
                            cubmat[i,j,k//nmo]['value'] = float(line[k%6-1])
        else:   # read real-space function
            if imo != 0:
                Warning('no orbital in cube file but imo != 0')
            # read grid data
            for i in range(int(n1)):
                for j in range(int(n2)):
                    iline = 1
                    line = f.readline().split()
                    for k in range(1, int(n3)+1):
                        if k // 6 >= iline:
                            line = f.readline().split()
                            iline += 1
                        cubmat[i,j,k-1]['value'] = float(line[k%6-1])
                        #line[-1] is same as line[5]
        
        # assign coordinate information to each grid point
        for i in range(int(n1)):
            for j in range(int(n2)):
                for k in range(int(n3)):
                    cubmat[i, j, k]['x'] = orgx + (i) * v1x + \
                        (j) * v2x + (k) * v3x
                    cubmat[i, j, k]['y'] = orgy + (i) * v1y + \
                        (j) * v2y + (k) * v3y
                    cubmat[i, j, k]['z'] = orgz + (i) * v1z + \
                        (j) * v2z + (k) * v3z
            
    return atoms, cubmat, nmo

def cmplx_amp_plot(title, window, real, img, x, y, z):
    '''
    plot amplitude contour of scalar complex array.
    !!! precise and controllable isosurface.
    !!! for MO, alpha and beta components should be plotted separately.
    ----------
    title : TYPE string
    ---title of figure.
    window : TYPE integer
    ---window number of the graphic.
    real : TYPE np.array
    ---real part of function.
    img : TYPE np.array
    ---imaginary part of function.
    x : TYPE np.array
    ---x grid points.
    y : TYPE np.array
    ---y grid points.
    z : TYPE np.array
    ---z grid points.
    -------
    Returns None
    '''
    fig = mlab.figure(window, size=(600,600), fgcolor=(0,0,0), \
        bgcolor=(1,1,1))
    fig_try = mlab.figure(675)
    contours = 0.1
    # plot isosurface of amplitude
    # amplitude of MO not amplitude of linear combine of AO
    amplitude = np.sqrt(np.add(np.square(real), np.square(img)))
    try:
        mlab.contour3d(x, y, z, amplitude, figure=fig_try, contours=[contours])
    except Exception as e:
        print(f'Error: {e}')
        if str(e).find('Contour instance must be') >= 0:
            contours = float(str(e)[str(e).rfind('<=')+3 : \
                str(e).rfind('<=')+14]) * .5 
                # factor of 0.5 is usually appropriate
            contour = mlab.contour3d(x, y, z, amplitude, figure=fig, \
                contours=[contours], opacity=.2, colormap='Greys')
            print('------')
            print('ignore previous error')
        else:
            exit()
    else:
        contour = mlab.contour3d(x, y, z, amplitude, figure=fig, \
            contours=[contours], opacity=.2, colormap='Greys')
    finally:
        mlab.title(f'{title}(amplitude), isovalue={contours}')
        contour.actor.property.line_width = 1
        contour.actor.property.edge_visibility = True
    mlab.close(675)

def cmplx_pha_plot(title, window, real, img, x, y, z):
    '''
    plot amplitude(contour) and phase(color) of scalar complex array.
    !!! can't specify the value of amplitude isovalue,
    shape of MO is only for phase distribution.
    !!! for MO, alpha and beta components should be plotted separately.
    ----------
    title : TYPE string
    ---title of figure.
    window : TYPE integer
    ---window number of the graphic.
    real : TYPE np.array
    ---real part of function.
    img : TYPE np.array
    ---imaginary part of function.
    x : TYPE np.array
    ---x grid points.
    y : TYPE np.array
    ---y grid points.
    z : TYPE np.array
    ---z grid points.
    -------
    Returns None
    '''
    fig = mlab.figure(window, size=(600,600), fgcolor=(0,0,0), \
        bgcolor=(1,1,1))
    # isosurface of amplitude
    # amplitude of MO are not linear combine of amplitude of AO
    amplitude = np.sqrt(np.add(np.square(real), np.square(img)))
    src = mlab.pipeline.scalar_field(amplitude)
    # colour mapping of phase
    # use arctan2 to describe complex plane completly
    phase = np.arctan2(real,img)
    src.image_data.point_data.add_array(phase.T.ravel())
    src.image_data.point_data.get_array(1).name = 'phase'
    src.update()
    
    # contour of amplitude
    src2 = mlab.pipeline.set_active_attribute(src, point_scalars='scalar')
    contour = mlab.pipeline.contour(src2, figure=fig)
    # contour of phase based on contour of amplitude
    contour2 = mlab.pipeline.set_active_attribute(contour,
        point_scalars='phase', figure=fig)
    mlab.pipeline.surface(contour2, colormap='hsv', opacity=.6, \
        vmax=np.pi, vmin=-np.pi)
    mlab.colorbar(title='phase', orientation='vertical', nb_labels=5)
    mlab.title(f'{title}(phase)')

def hex2rgb(hexcolor):
    '''
    colour format conversion from Hex to RGB.
    ----------
    hexcolor : TYPE string
    ---Hex format of colour.
    -------
    Returns turple*1
    '''
    hexcolor = int(hexcolor, base=16) \
        if isinstance(hexcolor, str) else hexcolor
    rgb = ((hexcolor >> 16) & 0xff, (hexcolor >> 8) & 0xff, hexcolor & 0xff)
    return rgb

def atoms_plot(window, atom_list, atomx, atomy, atomz):
    '''
    plot spherical model of atoms.
    ----------
    window : TYPE integer
    ---window number of the graphic.
    atom_list : TYPE np.array
    ---atomic number of each atom.
    atomx : TYPE np.array
    ---x coordinate of each atom.
    atomy : TYPE np.array
    ---y coordinate of each atom.
    atomz : TYPE np.array
    ---z coordinate of each atom.
    -------
    Returns None
    '''
    # load standard data
    json_address = environ.get('TRESC') + '\scripts\PubChemElements_all.json'
    with open(json_address,'r',encoding='utf-8') as f:
        std_dat = jsload(f)
    # get information of each atom
    atom_color=[]
    atom_symbol=[]
    atom_radius=[]
    num_atom = 0
    for i in atom_list:
        num_atom += 1
        if i <= 0 or i > 118:
            exit('atomic number cannot be recognised.')
        for j in std_dat["Table"]["Row"]:
            if i == int(j['Cell'][0]):
                atom_symbol.append(j["Cell"][1])
                atom_color.append(j["Cell"][4])
                atom_radius.append(qcel.covalentradii.get(i)) # covalent (CSD)
                #atom_radius.append(float(j["Cell"][7])/150.) # Van der Waal
                break
    # plot atoms
    fig = mlab.figure(window)
    for i in range(num_atom):
        RGB = hex2rgb(atom_color[i])
        R = RGB[0] / 255.
        G = RGB[1] / 255.
        B = RGB[2] / 255.
        mlab.points3d(atomx[i], atomy[i], atomz[i], \
            scale_factor=atom_radius[i], resolution=20, \
                color=(R, G, B), figure=fig, opacity=1., scale_mode='none')
        #mlab.text3d(atomx[i], atomy[i], atomz[i], \
        #    atom_symbol[i].strip(), line_width=.7, \
        #        color=(0,0,0), figure=fig, scale=(.3,.3,.3))
    

def _2corb(real_cube_address, img_cube_address):
    '''
    2-component complex MO visualization.
    ----------
    real_cube_address : TYPE string
    ---address of .cube file contains real part of MO.
    img_cube_address : TYPE string
    ---address of .cube file contains imaginary part of MO.
    -------
    Returns None
    '''
    atoms, real_mat_alpha, nmo = read_cube(real_cube_address, 1)
    if nmo == 0:
        exit('no orbital info in cube file.')
    elif nmo != 2:
        exit('cube file should contain both alpha and beta orbitals.')
    # grid data
    x = real_mat_alpha['x']
    y = real_mat_alpha['y']
    z = real_mat_alpha['z']
    real_val_alpha = real_mat_alpha['value']
    atoms, real_mat_beta, nmo = read_cube(real_cube_address, 2)
    real_val_beta = real_mat_beta['value']
    atoms, img_mat_alpha, nmo = read_cube(img_cube_address, 1)
    img_val_alpha = img_mat_alpha['value']
    atoms, img_mat_beta, nmo = read_cube(img_cube_address, 2)
    img_val_beta = img_mat_beta['value']
    
    # plot the amplitude of orbital
    cmplx_amp_plot('alpha', 1, real_val_alpha, img_val_alpha, x, y, z)
    cmplx_amp_plot('beta', 2, real_val_beta, img_val_beta, x, y, z)
    # plot the phase of orbital
    cmplx_pha_plot('alpha', 3, real_val_alpha, img_val_alpha, x, y, z)
    cmplx_pha_plot('beta', 4, real_val_beta, img_val_beta, x, y, z)
    
    # atoms data
    atom_list = atoms['index']
    atomx = atoms['x']
    atomy = atoms['y']
    atomz = atoms['z']
    
    # plot atoms, phase plot can't represent atomic position accurately
    atoms_plot(1, atom_list, atomx, atomy, atomz)
    atoms_plot(2, atom_list, atomx, atomy, atomz)
    
    mlab.show()

def run():
    real_cube_address = r'.\real.cub'
    img_cube_address = r'.\img.cub'
    _2corb(real_cube_address, img_cube_address)

if __name__ == '__main__':
    run()
