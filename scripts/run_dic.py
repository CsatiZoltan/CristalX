# -*- coding: utf-8 -*-
from os import path
from os.path import join

import numpy as np
import matplotlib.pyplot as plt

from grains import HAS_TABLES
from grains.geometry import TriMesh
from grains.dic import DIC, plot_strain

show_plots = True  # set to False to prevent from showing the plots

# Directory where we read from and write to
script_dir = path.dirname(__file__)  # absolute path to the directory the script is in
data_dir = join(script_dir, 'data')

# Options (https://stackoverflow.com/a/18633695/4892892,
# https://tex.stackexchange.com/a/391078/119426)
params = {'text.usetex': True,
          'text.latex.preamble': [r'\usepackage{siunitx}'],  # just an example to include a package
          'font.size': 16}
plt.rcParams.update(params)
# plt.rcParams

# Load the experimental data
# dic_snapshot = join(data_dir, '1_dic_field.npy')
dic_snapshot = join(data_dir, '2_dic_field.npy')
if path.isfile(dic_snapshot):
    u_restricted, v_restricted, _ = np.load(dic_snapshot, allow_pickle=True)
else:
    # Try to load it from the raw data
    # dic_data = join(data_dir, 'test.hdf')  # 6.3 GB, not committed to the Github repository
    dic_data = join(data_dir, 'FULLTEST_fields.hdf')  # 35.6 GB, not committed to the Github
    # repository
    if path.isfile(dic_data):
        if HAS_TABLES:
            # Fetch the displacement field from the series of field measurements
            import tables
            h = tables.open_file(dic_data)
            measurement = h.root.res
            time = 350  # out of 400
            u = measurement[time, :, :, 0]
            v = measurement[time, :, :, 1]

            # Show the axial strain
            if show_plots:
                plt.matshow(np.gradient(u)[1], interpolation='none', vmin=0, vmax=0.2)
                cbar = plt.colorbar(orientation='horizontal', format='%.2f', aspect=100,
                                    label=r'$\varepsilon_{xx}$')
                cbar.ax.set_xlabel(r'$\varepsilon_{xx}$', fontsize=30)

            # Crop data which belongs to the testing device (the indices were obtained using:
            # `selected = plt.ginput(-1, timeout=-1))`
            # restricted_region = np.s_[time, 111:914, 5:1963, :]
            restricted_region = np.s_[time, 800:2175, 860:3000, :]
            u_component = np.s_[:, :, 0]
            v_component = np.s_[:, :, 1]
            u_restricted = measurement[restricted_region][u_component]
            v_restricted = measurement[restricted_region][v_component]

            # Save the displacement components
            displacement_field = (u_restricted, v_restricted, 'Time step: 350 out of 400, '
                                  'cropped: 111:914, 5:1963 out of dimensions 1024x2048.')
            # np.save(join(data_dir, '1_dic_field.npy'), displacement_field)
            np.save(join(data_dir, '2_dic_field.npy'), displacement_field)
        else:
            raise ImportError('The PyTables package is needed to load the dataset.')
    else:
        raise Exception('Measurement dataset {0} not found.'.format(dic_data))
dic = DIC(u_restricted, v_restricted)
if show_plots:
    # Show the axial strain on the cropped region
    infinitesimal_strain = dic.strain('infinitesimal')
    plot_strain(infinitesimal_strain[:, :, 0], colorbar=True, label=r'$\varepsilon_{xx}$')
    plt.show()

# Set the scale, i.e. the physical length corresponding to the image data
dic.set_transformation((0, 0), 45.18)
dic.plot_physicalgrid()
plt.show()

# Superimpose the meshed domain of the numerical problem and the DIC grid
# mesh_file = join(data_dir, '1_mesh_extended_scaled.npz')
# with np.load(mesh_file, allow_pickle=True) as mesh:
#     if np.alltrue([variable_name in {'nodes', 'elements', 'element_groups', 'node_groups'} for
#                    variable_name in mesh]):
#         nodes = mesh['nodes']
#         elements = mesh['elements']
#         element_groups = mesh['element_groups']
#         node_groups = mesh['node_groups']
#     else:
#         raise Exception('Mesh must contain variables "nodes", "elements", "element_groups" and '
#                         '"node_groups".')
# mesh = TriMesh(nodes, elements)
# ax = dic.plot_superimposedmesh(mesh)
# ax.set_xlim(17.25, 19.25)
# ax.set_ylim(-10.8, -9.4)
# plt.show()
