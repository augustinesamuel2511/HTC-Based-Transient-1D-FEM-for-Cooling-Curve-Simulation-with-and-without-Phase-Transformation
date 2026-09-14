import numpy as np
from scipy import sparse
from scipy.sparse.linalg import cg, spsolve
import os


def load_preprocess_data(output_folder):
    data = {}
    # Coordinates
    cord_file = os.path.join(output_folder, "cordM.txt").replace('\\', '/')
    data["cordM"] = np.loadtxt(cord_file)

    # Connectivity
    con_file = os.path.join(output_folder, "conM.txt").replace('\\', '/')
    data["conM"] = np.loadtxt(con_file, dtype=int)

    # Coordinate system (string)
    cs_file = os.path.join(output_folder, "coordinate_system.txt").replace('\\', '/')
    with open(cs_file, "r") as f:
        data["coordinate_system"] = f.read().strip()

    pt_file = os.path.join(output_folder, "PT_file.txt").replace('\\', '/')
    with open(pt_file, "r") as f:
        pt_flag = f.read().strip()

    if pt_flag == "withPT":
        # density
        density_file = os.path.join(output_folder, "density_PT.txt").replace('\\', '/')
        data["density"] = np.loadtxt(density_file, dtype=float)

        # thermal conductivity
        k_file = os.path.join(output_folder, "thermal_conductivity_PT.txt").replace('\\', '/')
        data["thermal_conductivity"] = np.loadtxt(k_file, dtype=float)

        # specific heat
        cp_file = os.path.join(output_folder, "specific_heat_PT.txt").replace('\\', '/')
        data["specific_heat"] = np.loadtxt(cp_file, dtype=float)

        # latent heat
        h_file = os.path.join(output_folder, "latent_heat_PT.txt").replace('\\', '/')
        data["latent_heat"] = np.loadtxt(h_file, dtype=float)

        # TTT diagram
        ttt_file = os.path.join(output_folder, "ttt_diagram.txt").replace('\\', '/')
        data["ttt"] = np.loadtxt(ttt_file, dtype=float)

        # critical temperature
        crit_file = os.path.join(output_folder, "critical_temperatures.txt").replace('\\', '/')
        data["crit_temp"] = np.loadtxt(crit_file, dtype=float)

        # Initial Austenite
        in_aus_file = os.path.join(output_folder, "initial_austenite.txt").replace('\\', '/')
        data["init_aus"] = np.loadtxt(in_aus_file, dtype=float)

        # Equilibrium Ferrite
        femax_file = os.path.join(output_folder, "equilibrium_ferrite.txt").replace('\\', '/')
        data["Femax"] = np.loadtxt(femax_file, dtype=float)

    elif pt_flag == "noPT":
        mat_file = os.path.join(output_folder, "material_properties.txt").replace('\\', '/')
        data["material_properties"] = np.loadtxt(mat_file)
    return data


def mat_prop(matfile, temp):
    matfile = np.asarray(matfile)
    temp = np.asarray(temp).reshape(-1)

    # Define polynomials
    rho = lambda T: matfile[0, 0] * T ** 3 + matfile[0, 1] * T ** 2 + matfile[0, 2] * T + matfile[0, 3]
    K = lambda T: matfile[1, 0] * T ** 3 + matfile[1, 1] * T ** 2 + matfile[1, 2] * T + matfile[1, 3]
    Cp = lambda T: matfile[2, 0] * T ** 3 + matfile[2, 1] * T ** 2 + matfile[2, 2] * T + matfile[2, 3]

    k = K(temp)
    c = rho(temp) * Cp(temp)

    return c, k


def shapefn1D(gp):
    N = np.zeros((len(gp), 2))
    for i, xi in enumerate(gp):
        N[i, 0] = (1 - xi) / 2
        N[i, 1] = (1 + xi) / 2
    return N


def CKF(conM, cordM, lam, c, htc, coord_sys, Tamb):
    nNode_elem = 2
    dof_node = 1
    dof_elem = nNode_elem * dof_node
    nElems = conM.shape[0]
    nnodes = cordM.shape[0]
    nDOFs = nnodes * dof_node
    rmax = np.max(cordM) / 1000.0
    wt = 1  # for cartesian
    # Initialize global matrices
    Kg = np.zeros((nDOFs, nDOFs))
    Cg = np.zeros((nDOFs, nDOFs))
    Fg = np.zeros(nDOFs)

    # Gauss points and weights (3-point integration)
    gp1D = np.array([-np.sqrt(3 / 5), 0.0, np.sqrt(3 / 5)])
    w1D = np.array([5 / 9, 8 / 9, 5 / 9])

    # Shape functions and derivatives
    N1D = shapefn1D(gp1D)
    dN1D = np.array([[-0.5, 0.5]])  # constant for linear element

    # Loop over elements
    for k in range(nElems):
        Ke = np.zeros((dof_elem, dof_elem))
        Ce = np.zeros((dof_elem, dof_elem))
        Fe = np.zeros(dof_elem)

        # Coordinates of nodes of the element
        cordXk = cordM[conM[k, :]] / 1000.0

        for p in range(len(gp1D)):
            Xgp = N1D[p, :] @ cordXk
            J = (cordXk[1] - cordXk[0]) / 2.0
            if coord_sys == "Cylindrical":
                wt = Xgp
            for i in range(dof_elem):
                for j in range(dof_elem):
                    Ce[i, j] += c[conM[k, i]] * wt * N1D[p, i] * N1D[p, j] * J * w1D[p]
                    Ke[i, j] += lam[conM[k, i]] * wt * dN1D[0, i] * dN1D[0, j] * (1 / J) * w1D[p]

        # Assembly into global matrices
        for i in range(dof_elem):
            ig = conM[k, i]
            Fg[ig] += Fe[i]
            for j in range(dof_elem):
                jg = conM[k, j]
                Cg[ig, jg] += Ce[i, j]
                Kg[ig, jg] += Ke[i, j]

    # Apply unknown heat flux boundary condition
    if coord_sys == "Cartesian":
        Kg[-1, -1] += htc
        Fg[-1] += htc * Tamb

    elif coord_sys == "Cylindrical":
        Kg[-1, -1] += rmax * htc
        Fg[-1] += rmax * htc * Tamb

    Kg = sparse.csr_matrix(Kg)
    Cg = sparse.csr_matrix(Cg)
    Fg = np.asarray(Fg).ravel()

    return Cg, Kg, Fg


def FEM1D(Tin, Tamb, htc, dt, output_folder):
    meshdata = load_preprocess_data(output_folder)
    conM = meshdata['conM']
    cordM = meshdata['cordM']
    coord_sys = meshdata['coordinate_system']
    matfile = meshdata['material_properties']
    theta = 0.6667
    nNodes = cordM.shape[0]

    T = np.zeros((nNodes, 1))
    T[:, 0] = Tin

    for ts in range(1, 20000000):
        T_temp = T[:, ts - 1].copy()
        ep = 100.0

        # Compute initial C, K, F
        c_arr, k_arr = mat_prop(matfile, T_temp)
        C, K, F = CKF(conM, cordM, k_arr, c_arr, htc, coord_sys, Tamb)
        B_n = (C - dt * (1 - theta) * K) @ T_temp + dt * (1 - theta) * F

        ss = 1
        while ep > 1e-4:
            ss += 1
            c_arr, k_arr = mat_prop(matfile, T_temp)
            _, K_temp, F_temp = CKF(conM, cordM, k_arr, c_arr, htc, coord_sys, Tamb)

            A = (dt * theta) * K_temp + C
            B = B_n + (dt * theta) * F_temp

            T_new = spsolve(A, B)

            # Prevent temperature from falling below ambient
            T_new = np.maximum(T_new, Tamb)

            ep = np.linalg.norm(T_new - T_temp.ravel()) / (np.linalg.norm(T_temp) + 1e-16)
            T_temp = T_new

            if ss > 200:
                break

        T = np.column_stack((T, T_temp))
        print(f"\ntime step = {ts}")
        if T_temp[0] <= Tamb + 0.01:
            T_temp[:] = Tamb
            T[:, ts] = T_temp
            break
    # --- Save temperature matrix T as TXT with tab delimiter ---
    T_path = os.path.join(output_folder, 'temperature_matrix.txt').replace('\\', '/')

    np.savetxt(T_path, np.round(T, 2), delimiter='\t')
    print(f"Temperature matrix saved to {T_path}")
