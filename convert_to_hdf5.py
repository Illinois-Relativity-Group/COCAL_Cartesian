import h5py
import numpy as np

data = np.loadtxt("data.txt", skiprows=1)
x, y, z, rho_ns, rho_disk  = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:,4]

x_unique = np.unique(x)
y_unique = np.unique(y)
z_unique = np.unique(z)

##print("Max density is: " + str(np.max(density)))
print("Max rho_disk:", np.max(rho_disk), "Max rho_ns:", np.max(rho_ns))
grid_shape = (len(x_unique), len(y_unique), len(z_unique))

##density_grid = density.reshape(grid_shape)
density_grid_disk = rho_disk.reshape(grid_shape)
density_grid_ns = rho_ns.reshape(grid_shape)

output_path = 'plot_data.h5'
with h5py.File(output_path, 'w') as hdf:
    hdf.attrs['grid_type'] = 'structured'
    hdf.attrs['x_range'] = (x_unique[0], x_unique[-1])
    hdf.attrs['y_range'] = (y_unique[0], y_unique[-1])
    hdf.attrs['z_range'] = (z_unique[0], z_unique[-1])
    hdf.attrs['refinement_levels'] = 1  # Single level, no refinement
    hdf.create_dataset('rho_disk', data=density_grid_disk, dtype='f4')
    hdf.create_dataset('rho_ns', data=density_grid_ns, dtype='f4')
    ##hdf.create_dataset('density', data=density_grid, dtype='f4')
    hdf.create_dataset('x', data=x_unique, dtype='f4')
    hdf.create_dataset('y', data=y_unique, dtype='f4')
    hdf.create_dataset('z', data=z_unique, dtype='f4')

print("File saved as: " + output_path)

with h5py.File(output_path, 'r') as hdf:
    print("File contents:")
    for key in hdf.keys():
        print(f"{key}: {hdf[key].shape}, dtype={hdf[key].dtype}")
