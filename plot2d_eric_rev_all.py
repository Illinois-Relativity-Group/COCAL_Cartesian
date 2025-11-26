import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Color maps for NS and disk
hex_list_ns = ['#bdbca2', '#bdbca2','#5500a4','#7803fb','#9309dd','#aa1657','#be2c00','#d04c00','#e17800','#f0b300','#ffff00']
hex_list_disk = ['#bdbca2', '#bdbca2', '#a8f79d', '#8dfead', '#72febc', '#57f7c9', '#3cead5', '#22d5e0', '#07bcea', '#149df1', '#2f79f7']
float_list=[0, 0.001, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]

def hex_to_rgb(value):
    value = value.strip("#")
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def get_cmap(hex_list, float_list=None):
    rgb_list = [[v/256. for v in hex_to_rgb(h)] for h in hex_list]
    cdict = {}
    for num, col in enumerate(['red', 'green', 'blue']):
        cdict[col] = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] 
                      for i in range(len(float_list))]
    return mcolors.LinearSegmentedColormap('custom', segmentdata=cdict, N=256)

cmap_ns = get_cmap(hex_list_ns, float_list)    
cmap_disk = get_cmap(hex_list_disk, float_list) 

print("Loading data...")
data = np.loadtxt('data.txt', skiprows=1)
x, y, z = data[:, 0], data[:, 1], data[:, 2]
rho_ns, rho_disk = data[:, 3], data[:, 4]

tolerance = 1
thresh = 0.001

# XY slice (z=0)
print("Creating XY slice...")
plt.figure(figsize=(14, 5))
mask_xy = np.abs(z - 0) < tolerance
xs_xy, ys_xy = x[mask_xy], y[mask_xy]
rho_ns_xy = np.where(rho_ns[mask_xy]/np.max(rho_ns) > thresh, rho_ns[mask_xy]/np.max(rho_ns), np.nan)
rho_disk_xy = np.where(rho_disk[mask_xy]/np.max(rho_disk) > thresh, rho_disk[mask_xy]/np.max(rho_disk), np.nan)

sc1 = plt.scatter(xs_xy, ys_xy, c=rho_ns_xy, cmap=cmap_ns, s=3, vmin=0, vmax=1, marker='s')
sc2 = plt.scatter(xs_xy, ys_xy, c=rho_disk_xy, cmap=cmap_disk, s=3, vmin=0, vmax=1, marker='s')
plt.xlabel('x', fontsize=12)
plt.ylabel('y', fontsize=12)
plt.title('XY slice (z=0)', fontsize=12)
plt.gca().set_aspect('equal')
plt.gca().set_facecolor('#bdbca2')
cbar1 = plt.colorbar(sc1, pad=0.02, aspect=30)
cbar1.set_label('NS density (normalized)', fontsize=10)
cbar2 = plt.colorbar(sc2, pad=0.12, aspect=30)
cbar2.set_label('Disk density (normalized)', fontsize=10)
plt.tight_layout()
plt.savefig('nst_slice_xy.png', dpi=300, bbox_inches='tight')
print("Saved nst_slice_xy.png")
plt.close()

# XZ slice (y=0)
print("Creating XZ slice...")
plt.figure(figsize=(14, 5))
mask_xz = np.abs(y - 0) < tolerance
xs_xz, zs_xz = x[mask_xz], z[mask_xz]
rho_ns_xz = np.where(rho_ns[mask_xz]/np.max(rho_ns) > thresh, rho_ns[mask_xz]/np.max(rho_ns), np.nan)
rho_disk_xz = np.where(rho_disk[mask_xz]/np.max(rho_disk) > thresh, rho_disk[mask_xz]/np.max(rho_disk), np.nan)

sc1 = plt.scatter(xs_xz, zs_xz, c=rho_ns_xz, cmap=cmap_ns, s=3, vmin=0, vmax=1, marker='s')
sc2 = plt.scatter(xs_xz, zs_xz, c=rho_disk_xz, cmap=cmap_disk, s=3, vmin=0, vmax=1, marker='s')
plt.xlabel('x', fontsize=12)
plt.ylabel('z', fontsize=12)
plt.title('XZ slice (y=0)', fontsize=12)
plt.gca().set_aspect('equal')
plt.gca().set_facecolor('#bdbca2')
cbar1 = plt.colorbar(sc1, pad=0.02, aspect=30)
cbar1.set_label('NS density (normalized)', fontsize=10)
cbar2 = plt.colorbar(sc2, pad=0.12, aspect=30)
cbar2.set_label('Disk density (normalized)', fontsize=10)
plt.tight_layout()
plt.savefig('nst_slice_xz.png', dpi=300, bbox_inches='tight')
print("Saved nst_slice_xz.png")
plt.close()

# YZ slice (x=0)
print("Creating YZ slice...")
plt.figure(figsize=(14, 5))
mask_yz = np.abs(x - 0) < tolerance
ys_yz, zs_yz = y[mask_yz], z[mask_yz]
rho_ns_yz = np.where(rho_ns[mask_yz]/np.max(rho_ns) > thresh, rho_ns[mask_yz]/np.max(rho_ns), np.nan)
rho_disk_yz = np.where(rho_disk[mask_yz]/np.max(rho_disk) > thresh, rho_disk[mask_yz]/np.max(rho_disk), np.nan)

sc1 = plt.scatter(ys_yz, zs_yz, c=rho_ns_yz, cmap=cmap_ns, s=3, vmin=0, vmax=1, marker='s')
sc2 = plt.scatter(ys_yz, zs_yz, c=rho_disk_yz, cmap=cmap_disk, s=3, vmin=0, vmax=1, marker='s')
plt.xlabel('y', fontsize=12)
plt.ylabel('z', fontsize=12)
plt.title('YZ slice (x=0)', fontsize=12)
plt.gca().set_aspect('equal')
plt.gca().set_facecolor('#bdbca2')
cbar1 = plt.colorbar(sc1, pad=0.02, aspect=30)
cbar1.set_label('NS density (normalized)', fontsize=10)
cbar2 = plt.colorbar(sc2, pad=0.12, aspect=30)
cbar2.set_label('Disk density (normalized)', fontsize=10)
plt.tight_layout()
plt.savefig('nst_slice_yz.png', dpi=300, bbox_inches='tight')
print("Saved nst_slice_yz.png")
plt.close()

print("All done!")