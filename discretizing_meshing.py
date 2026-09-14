import numpy as np
import os


def discretizing_meshing(l, nx, output_folder):
    lnodes = 2  # number of nodes per element (local nodes)
    npx = nx + 1  # total number of points in x

    # Discretizing length
    coords = np.linspace(0, l, npx)

    cordM = coords.reshape(-1, 1)  # coordinate matrix as column vector
    gnodes = np.arange(npx)  # global nodes

    # Initialize connectivity matrix
    conM = np.zeros((nx, lnodes), dtype=int)
    conM[:, 0] = gnodes[:-1]
    conM[:, 1] = gnodes[1:]

    # --- Save to separate files if output_folder is provided ---
    if output_folder:
        os.makedirs(output_folder, exist_ok=True)

        np.savetxt(os.path.join(output_folder, "cordM.txt").replace('\\', '/'), cordM, fmt="%.6f")
        np.savetxt(os.path.join(output_folder, "conM.txt").replace('\\', '/'), conM, fmt="%d")
