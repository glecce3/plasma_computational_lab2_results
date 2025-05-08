import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LogNorm, Normalize, LinearSegmentedColormap
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.constants import micron, c, pi, centi, femto, e, epsilon_0, m_e, m_p
import happi
import os
import cv2
# Use Agg backend for matplotlib
matplotlib.use('Agg')

# Constants
lambda_SI = 0.8 * micron
omega_SI = 2.0 * pi * c / lambda_SI
mc2 = 0.510998950e6
n_crit = (m_e * epsilon_0 * (2 * pi * c)**2) / ((e * lambda_SI)**2)  # Critical plasma density
my_dpi = 300  # DPI for figures

# Open Smilei simulation data
smilei_dir = '.'
s = happi.Open(smilei_dir)

# Simulation parameters
um_s = 1 / s.namelist.um  # Conversion to micron
fs_s = 1 / omega_SI * 1e15  # Conversion to fs
Lx = s.namelist.Lx * um_s  # Simulation domain size
Ly = s.namelist.Ly * um_s
unitsim = ["fs", "um", "MeV", "J/m"]  # Units for simulation output
min_density_plot_density = 0.01 # units of n_c
max_density_plot_density = 1 # units of n_c
Z_ion = 1

# Get colormap
ncolors = 256
color_array = plt.get_cmap('seismic')(range(ncolors))
color_array[:, -1] = np.concatenate(
    (np.ones(int(14 * ncolors / 32)), np.zeros(int(4 * ncolors / 32)), np.ones(int(14 * ncolors / 32))))
map_object = LinearSegmentedColormap.from_list(name='diverg', colors=color_array)
matplotlib.colormaps.register(cmap=map_object)

# Iterate over timesteps
data = s.Field(0, 'Ex')
steps = data.getTimesteps()

for n, ts in enumerate(steps):
    fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(10, 5), dpi=my_dpi, sharex=True, sharey=True)

    # Electron density and Bzf
    ne = s.Field(0, "Rho_ele", units=unitsim)
    ni = s.Field(0, "Rho_ion", units=unitsim)
    time = ne.getTimes()[n]
    x = ne.getAxis('x')
    y = ne.getAxis('y')
    rhoele = -ne.getData(ts)[0]
    rhoion = ni.getData(ts)[0]
    bz = s.Field(0, "Bz", units=unitsim)
    Bz = bz.getData(timestep=ts)[0]

    # Plot electron density and Bz
    ax[0].set_title('Ele+Bz, t = %.2f fs' % time)
    ax[1].set_title('Ions, t = %.2f fs' % time)

    im = ax[0].imshow(np.flip(rhoele, axis=1).T, cmap='Greys', norm=LogNorm(vmin=min_density_plot_density, vmax=max_density_plot_density),
                      extent=[x.min(), x.max(), y.min(), y.max()], interpolation='bicubic')
    divider = make_axes_locatable(ax[0])
    cax = divider.append_axes('bottom', size='3%', pad=0)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'n$_e$ [n$_c$]')
    im = ax[0].imshow(np.flip(Bz, axis=1).T, cmap='diverg', extent=[x.min(), x.max(), y.min(), y.max()],
                      vmin=-10, vmax=10, interpolation='bicubic')
    divider = make_axes_locatable(ax[0])
    cax = divider.append_axes('bottom', size='3%', pad=1.1)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=r'B$_z$ [T]')

    # Plot ion density
    im = ax[1].imshow(np.flip(rhoion, axis=1).T, cmap='Blues', norm=LogNorm(vmin=min_density_plot_density/Z_ion, vmax=max_density_plot_density/Z_ion),
                      extent=[x.min(), x.max(), y.min(), y.max()], interpolation='bicubic')
    divider = make_axes_locatable(ax[1])
    cax = divider.append_axes('bottom', size='3%', pad=1.1)
    plt.colorbar(im, cax=cax, orientation='horizontal', label=r'n$_i$ [n$_c$]')

    # Set labels
    for a in ax.reshape(-1):
        a.set_xlabel(r'x [$\mu$m]')
        a.set_ylabel(r'y [$\mu$m]')


    # Specify your folder name (can also be an absolute path if needed)
    output_folder = os.path.expanduser('~/Desktop/output_maps')
    
    # Create the folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Compose the full file path
    image_file_name = os.path.join(output_folder, 'maps_%04d.png' % ts)
    
    # Save the figure
    plt.tight_layout()
    plt.savefig(image_file_name, dpi=my_dpi)
    plt.close()
    
    

# Folder containing your images
image_folder = os.path.expanduser('~/Desktop/output_maps')
video_name = os.path.join(image_folder, 'output_video.mp4')
fps = 6  # frames per second

# Get list of image files sorted (make sure naming allows correct sorting)
images = sorted([img for img in os.listdir(image_folder) if img.endswith(".png")])

# Read first image to get frame size
frame_path = os.path.join(image_folder, images[0])
frame = cv2.imread(frame_path)
height, width, layers = frame.shape

# Define the codec and create VideoWriter object
fourcc = cv2.VideoWriter_fourcc(*'mp4v')  # for .mp4 output
video = cv2.VideoWriter(video_name, fourcc, fps, (width, height))

# Add all images to video
for image in images:
    img_path = os.path.join(image_folder, image)
    frame = cv2.imread(img_path)
    video.write(frame)

video.release()
print("Video saved as", video_name)
