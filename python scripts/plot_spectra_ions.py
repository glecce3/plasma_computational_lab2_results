import happi
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
from scipy.constants import micron, c, pi, e, epsilon_0, m_e
import cv2

# Use Agg backend for matplotlib
matplotlib.use('Agg')

# Constants and quantities
wavelength = 0.8e-6  # Wavelength in meters
omega = 2 * pi * c / wavelength  # Laser frequency
nc = (m_e * epsilon_0 * (2 * pi * c)**2) / ((e * wavelength)**2)  # Critical plasma density

unitsim = ['MeV', 's']  # Units for simulation data

# Open Smilei simulation
S = happi.Open('.', reference_angular_frequency_SI=omega)

# Particle binning diagnostic
Diag = S.ParticleBinning(0, units=unitsim)
timesteps = Diag.getAvailableTimesteps()
times = Diag.getTimes()
E = Diag.getAxis('ekin')  # Kinetic energy axis

# Iterate over timesteps
for n, ts in enumerate(timesteps[::10]):
    # Calculate electron energy spectrum
    Wele = Diag.getData(ts)[0] * nc * (S.namelist.Lx * S.namelist.Ly * (c / omega)**2) / (
                m_e * c**2 / e * 1e-6)  # Electron energy spectrum
    fig, ax = plt.subplots(1, 1, figsize=(4, 3.5), dpi=300)
    time = times[n]  # Current time
    ax.plot(E, Wele, label=r'Ions', color='blue')  # Plot electron spectrum
    ax.set_xlabel(r'$E_p$ [MeV]')  # Set x-axis label
    leg = ax.legend(loc='best', ncol=1, shadow=False, fancybox=False, frameon=False)  # Add legend
    ax.set_yscale(r'log')  # Set y-axis scale to logarithmic
    ax.set_ylabel(r'$dN/dE$ [arb. units]')  # Set y-axis label
    ax.set_title(f't = {time * 1e16:.2f} fs')   
    plt.tight_layout()
    plt.tight_layout()
    output_folder = os.path.expanduser('~/Desktop/output_spectra/ions')
    os.makedirs(output_folder, exist_ok=True)
    file_name = 'spectra_%04d.png' % ts
    full_path = os.path.join(output_folder, file_name)
# Save and close figure
    fig.savefig(full_path)
    plt.close()


# Folder containing your images
image_folder = os.path.expanduser('~/Desktop/output_spectra/ions')
video_name = os.path.join(image_folder, 'output_video_spectra_ions.mp4')
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

Wele = Diag.getData(2130)[0] * nc * (S.namelist.Lx * S.namelist.Ly * (c / omega)**2) / (
            m_e * c**2 / e * 1e-6)  # Electron energy spectrum; 2130 => 301.44fs;
# Supponiamo che E e Wele siano array 1D della stessa lunghezza
data = np.column_stack((E, Wele))  # Combina E e Wions in due colonne

# Salva in CSV con intestazione
np.savetxt('spectra_ions.csv', data, delimiter=',', header='Eion,Wion', comments='')
