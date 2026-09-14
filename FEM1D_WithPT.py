import numpy as np
import pandas as pd
from scipy.interpolate import pchip_interpolate
from scipy import sparse
from scipy.sparse.linalg import spsolve
from FEM1D_NoPT import load_preprocess_data, shapefn1D
import os


def SchielSum(S, A, T, dt):
    # t = interp1(A(:,1), A(:,2), T, 'pchip', 'extrap')
    t = pchip_interpolate(A[:, 0], A[:, 1], T)

    # s = ts / t
    s = S + (dt / t)

    return s


def b_n_calc(A, T):
    x1 = 0.1 / 100
    x2 = 99.9 / 100

    # Column mapping similar to MATLAB: A(:,1), A(:,2), A(:,3)
    t1 = pchip_interpolate(A[:, 0], A[:, 1], T)
    t2 = pchip_interpolate(A[:, 0], A[:, 2], T)

    n = np.log(np.log(1 - x1) / np.log(1 - x2)) / np.log(t1 / t2)
    b = -np.log(1 - x1) / (t1 ** n)

    return b, n


def fraction_transformed(T, Xk, A, dt):
    b, n = b_n_calc(A, T)
    eps = np.finfo(float).eps
    Xk_clamped = np.clip(Xk, 0.0, 1.0 - eps)

    # tj = ( -log(1-Xk) / b )^(1/n)
    tj = (-np.log(1.0 - Xk_clamped) / b) ** (1.0 / n)

    t = tj + dt

    X = 1.0 - np.exp(-b * (t ** n))

    # clamp to [0,1]
    X = np.clip(X, 0.0, 1.0)

    return X


def Transformation(T, XM, XB, XP, XF, Sf, Sp, Sb, CritTemp,
                   dt, BTnb, FTnb, PTnb, in_aus, femax):
    T = np.asarray(T).ravel()
    nn = len(T)

    # Separate Schiel sums
    sf = Sf.copy()
    sp = Sp.copy()
    sb = Sb.copy()

    Ae3, Ae1, Bs, Ms = CritTemp[0], CritTemp[1], CritTemp[2], CritTemp[3]

    Xf = XF.copy()
    Xp = XP.copy()
    Xb = XB.copy()
    Xm = XM.copy()

    # Austenite fraction
    XA = in_aus - XF - XP - XB

    for i in range(nn):

        # ---------------------------------------------------------
        # Martensite
        # ---------------------------------------------------------
        if T[i] < Ms:

            Xm[i] = XA[i] * (1 - np.exp(-0.011 * (Ms - T[i])))

            if XM[i] > Xm[i]:
                Xm[i] = XM[i]

        # ---------------------------------------------------------
        # Bainite
        # ---------------------------------------------------------
        elif Bs > T[i] >= Ms:

            sb[i] = SchielSum(Sb[i], BTnb, T[i], dt)

            if sb[i] >= 1:

                if XB[i] == 0:
                    XB[i] = 0.001

                den = max((in_aus - Xf[i] - Xp[i]) * (in_aus - (XB[i] + Xf[i] + Xp[i])), 1e-12)

                B = XB[i] / den

                Xb_temp = fraction_transformed(T[i], B, BTnb, dt)

                Xb[i] = Xb_temp * (in_aus - Xf[i] - Xp[i]) * (in_aus - (XB[i] + Xf[i] + Xp[i]))

        # ---------------------------------------------------------
        # Pearlite
        # ---------------------------------------------------------
        elif Ae1 >= T[i] >= Bs:

            sp[i] = SchielSum(Sp[i], PTnb, T[i], dt)

            if sp[i] >= 1:

                if XP[i] == 0:
                    XP[i] = 0.001

                den = max((in_aus - Xf[i]) * (in_aus - (XP[i] + Xf[i])), 1e-12)

                P = XP[i] / den

                Xp_temp = fraction_transformed(T[i], P, PTnb, dt)

                Xp[i] = Xp_temp * (in_aus - Xf[i]) * (in_aus - (XP[i] + Xf[i]))

        # ---------------------------------------------------------
        # Ferrite
        # ---------------------------------------------------------
        elif Ae3 != 0 and Ae3 >= T[i] > Ae1:

            sf[i] = SchielSum(Sf[i], FTnb, T[i], dt)

            if sf[i] >= 1:

                if XF[i] == 0:
                    XF[i] = 0.001

                Femax = femax * ((Ae3 - T[i]) / (Ae3 - Ae1))

                F = XF[i] / (Femax * (in_aus - XF[i]))

                Xf_temp = fraction_transformed(T[i], F, FTnb, dt)

                Xf[i] = Xf_temp * Femax * (in_aus - XF[i])

    return Xm, Xb, Xp, Xf, sf, sp, sb


def ThermalPropEnthalpy(T_temp, XM, XB, XP, XF,
                        Sf, Sp, Sb, dt, data):
    # ============================================================
    # Material property tables
    # ============================================================

    Tprop_rho = data['density'][:, 0]
    rho_prop = data['density'][:, 1:]

    Tprop_k = data['thermal_conductivity'][:, 0]
    kprop = data['thermal_conductivity'][:, 1:]

    Tprop_cp = data['specific_heat'][:, 0]
    cp_prop = data['specific_heat'][:, 1:]

    Tprop_h = data['latent_heat'][:, 0]
    Hprop = data['latent_heat'][:, 1:]

    crit_temp = data['crit_temp']

    # ============================================================
    # Distinguish Hypereutectoid / Hypoeutectoid steel
    #
    # crit_temp[0] == 0  --> Hypereutectoid
    # crit_temp[0] != 0  --> Hypoeutectoid
    # ============================================================

    if crit_temp[0] == 0:

        # ========================================================
        # HYPEREUTECTOID STEEL
        #
        # Phases:
        # Austenite
        # Pearlite
        # Bainite
        # Martensite
        #
        # Ferrite transformation is disabled.
        # ========================================================

        ftnb = None

        # -------------------------------
        # TTT data
        # -------------------------------

        # Pearlite: columns 0, 1, 2
        ptnb = np.column_stack((
            data['ttt'][:, 0],
            data['ttt'][:, 1],
            data['ttt'][:, 2]
        ))

        # Bainite: columns 0, 3, 4
        btnb = np.column_stack((
            data['ttt'][:, 0],
            data['ttt'][:, 3],
            data['ttt'][:, 4]
        ))

        # -------------------------------
        # Property interpolation
        # -------------------------------

        def prop_estimate(Tprop, prop, T):

            # [Austenite, Pearlite, Bainite, Martensite]

            return np.column_stack((
                pchip_interpolate(Tprop, prop[:, 0], T),
                pchip_interpolate(Tprop, prop[:, 1], T),
                pchip_interpolate(Tprop, prop[:, 2], T),
                pchip_interpolate(Tprop, prop[:, 3], T)
            ))

        def latent_heat_estimate(Tprop, prop, T):

            # [Pearlite, Bainite, Martensite]

            return np.column_stack((
                pchip_interpolate(Tprop, prop[:, 0], T),
                pchip_interpolate(Tprop, prop[:, 1], T),
                pchip_interpolate(Tprop, prop[:, 2], T)
            ))

    else:

        # ========================================================
        # HYPOEUTECTOID STEEL
        #
        # Phases:
        # Austenite
        # Ferrite
        # Pearlite
        # Bainite
        # Martensite
        # ========================================================

        # -------------------------------
        # TTT data
        # -------------------------------

        # Ferrite: columns 0, 1, 2
        ftnb = np.column_stack((
            data['ttt'][:, 0],
            data['ttt'][:, 1],
            data['ttt'][:, 2]
        ))

        # Pearlite: columns 0, 2, 3
        ptnb = np.column_stack((
            data['ttt'][:, 0],
            data['ttt'][:, 2],
            data['ttt'][:, 3]
        ))

        # Bainite: columns 0, 4, 5
        btnb = np.column_stack((
            data['ttt'][:, 0],
            data['ttt'][:, 4],
            data['ttt'][:, 5]
        ))

        # -------------------------------
        # Property interpolation
        # -------------------------------

        def prop_estimate(Tprop, prop, T):

            # [Austenite, Ferrite, Pearlite, Bainite, Martensite]

            return np.column_stack((
                pchip_interpolate(Tprop, prop[:, 0], T),
                pchip_interpolate(Tprop, prop[:, 1], T),
                pchip_interpolate(Tprop, prop[:, 2], T),
                pchip_interpolate(Tprop, prop[:, 3], T),
                pchip_interpolate(Tprop, prop[:, 4], T)
            ))

        def latent_heat_estimate(Tprop, prop, T):

            # [Ferrite, Pearlite, Bainite, Martensite]

            return np.column_stack((
                pchip_interpolate(Tprop, prop[:, 0], T),
                pchip_interpolate(Tprop, prop[:, 1], T),
                pchip_interpolate(Tprop, prop[:, 2], T),
                pchip_interpolate(Tprop, prop[:, 3], T)
            ))

    # ============================================================
    # Initial austenite fraction and maximum ferrite fraction
    # ============================================================

    in_aus = data["init_aus"]
    femax = data["Femax"]

    # ============================================================
    # Phase transformation
    # ============================================================

    XM_new, XB_new, XP_new, XF_new, Sf_new, Sp_new, Sb_new = \
        Transformation(
            T_temp,
            XM,
            XB,
            XP,
            XF,
            Sf,
            Sp,
            Sb,
            crit_temp,
            dt,
            btnb,
            ftnb,
            ptnb,
            in_aus,
            femax
        )

    # ============================================================
    # Calculate phase fractions and latent heat
    # ============================================================

    if crit_temp[0] == 0:

        # ========================================================
        # HYPEREUTECTOID
        #
        # Phase order:
        # [Austenite, Pearlite, Bainite, Martensite]
        # ========================================================

        XA = in_aus - XP_new - XB_new - XM_new

        Xprop = np.column_stack((
            XA,
            XP_new,
            XB_new,
            XM_new
        ))

        # Change in phase fraction
        dX = np.column_stack((
            XP_new - XP,
            XB_new - XB,
            XM_new - XM
        ))

        # Density:
        # [Austenite, Pearlite, Bainite, Martensite]
        rho_all = prop_estimate(
            Tprop_rho,
            rho_prop,
            T_temp
        )

        # Remove austenite
        rho_interp = rho_all[:, 1:]

        # Latent heat:
        # [Pearlite, Bainite, Martensite]
        L_interp = np.abs(
            latent_heat_estimate(
                Tprop_h,
                Hprop,
                T_temp
            )
        )

    else:

        # ========================================================
        # HYPOEUTECTOID
        #
        # Phase order:
        # [Austenite, Ferrite, Pearlite, Bainite, Martensite]
        # ========================================================

        XA = in_aus - XF_new - XP_new - XB_new - XM_new

        Xprop = np.column_stack((
            XA,
            XF_new,
            XP_new,
            XB_new,
            XM_new
        ))

        # Change in phase fraction
        dX = np.column_stack((
            XF_new - XF,
            XP_new - XP,
            XB_new - XB,
            XM_new - XM
        ))

        # Density:
        # [Austenite, Ferrite, Pearlite, Bainite, Martensite]
        rho_all = prop_estimate(
            Tprop_rho,
            rho_prop,
            T_temp
        )

        # Remove austenite
        rho_interp = rho_all[:, 1:]

        # Latent heat:
        # [Ferrite, Pearlite, Bainite, Martensite]
        L_interp = np.abs(
            latent_heat_estimate(
                Tprop_h,
                Hprop,
                T_temp
            )
        )

    # ============================================================
    # Latent heat generation
    # ============================================================

    H = np.sum(
        rho_interp * L_interp * dX / dt,
        axis=1
    )

    # ============================================================
    # Thermal conductivity
    # ============================================================

    k_all = prop_estimate(
        Tprop_k,
        kprop,
        T_temp
    )

    k = np.sum(
        k_all * Xprop,
        axis=1
    )

    # ============================================================
    # Volumetric heat capacity
    # rho * cp
    # ============================================================

    rho_all = prop_estimate(
        Tprop_rho,
        rho_prop,
        T_temp
    )

    cp_all = prop_estimate(
        Tprop_cp,
        cp_prop,
        T_temp
    )

    rhoc = rho_all * cp_all

    c = np.sum(
        rhoc * Xprop,
        axis=1
    )

    # ============================================================
    # Return
    # ============================================================

    return (
        H,
        k,
        c,
        XM_new,
        XB_new,
        XP_new,
        XF_new,
        Sf_new,
        Sp_new,
        Sb_new
    )


def CKF_withpt(conM, cordM, H, lam, c, htc, coord_sys, Tamb):
    nNode_elem = 2
    dof_node = 1
    dof_elem = nNode_elem * dof_node
    wt = 1.0

    nElems = conM.shape[0]
    nnodes = cordM.shape[0]
    nDOFs = nnodes * dof_node

    rmax = np.max(cordM) / 1000

    # Global matrices
    Kg = np.zeros((nDOFs, nDOFs))
    Cg = np.zeros((nDOFs, nDOFs))
    Fg = np.zeros(nDOFs)

    # 1D 3-point Gauss rule
    gp = np.array([-np.sqrt(3 / 5), 0.0, np.sqrt(3 / 5)])  # Gauss points
    w = np.array([5 / 9, 8 / 9, 5 / 9])  # weights

    # Shape functions and derivatives
    N1D = shapefn1D(gp)
    dN1D = np.array([[-0.5, 0.5]])  # constant for linear element

    # --- Element Loop ---
    for k in range(nElems):
        Ke = np.zeros((dof_elem, dof_elem))
        Ce = np.zeros((dof_elem, dof_elem))
        Fe = np.zeros(dof_elem)

        # Extract nodal coordinates for this element
        # Coordinates of nodes of the element
        cordXk = cordM[conM[k, :]] / 1000.0

        for p in range(len(gp)):
            Xgp = N1D[p, :] @ cordXk
            J = (cordXk[1] - cordXk[0]) / 2.0
            if coord_sys == "Cylindrical":
                wt = Xgp

            for i in range(dof_elem):
                for j in range(dof_elem):
                    Ce[i, j] += c[conM[k, i]] * wt * N1D[p, i] * N1D[p, j] * J * w[p]
                    Ke[i, j] += lam[conM[k, i]] * wt * dN1D[0, i] * dN1D[0, j] * (1 / J) * w[p]

                Fe[i] += H[conM[k, i]] * wt * N1D[p, i] * J * w[p]

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
    Fg = np.asarray(Fg).ravel()
    return Cg, Kg, Fg


def FEM1D_withpt(Tin, Tamb, htc, dt, output_folder):
    data = load_preprocess_data(output_folder)
    conM = data['conM']
    cordM = data['cordM']
    coord_sys = data['coordinate_system']
    in_aus = data["init_aus"]
    theta = 0.6667
    nNodes = cordM.shape[0]

    # Temperature matrix
    T = np.zeros((nNodes, 1))
    T[:, 0] = Tin
    # Initialize phase fraction arrays
    XA = np.zeros((nNodes, 1))
    XA[:, 0] = in_aus
    XF = np.zeros((nNodes, 1))  # Ferrite
    XM = np.zeros((nNodes, 1))  # Martensite
    XB = np.zeros((nNodes, 1))  # Bainite
    XP = np.zeros((nNodes, 1))  # Pearlite
    Sf = np.zeros((nNodes, 1))  # Ferrite Schiel sum
    Sp = np.zeros((nNodes, 1))  # Pearlite Schiel sum
    Sb = np.zeros((nNodes, 1))  # Bainite Schiel sum

    # --- Time stepping ---
    for ts in range(1, 20000000):
        T_temp = T[:, ts - 1].copy()

        # Compute initial properties
        ep = 100
        H, k, c, xm, xb, xp, xf, Sf, Sp, Sb = ThermalPropEnthalpy(
            T_temp,
            XM[:, ts - 1],
            XB[:, ts - 1],
            XP[:, ts - 1],
            XF[:, ts - 1],
            Sf, Sp, Sb,
            dt,
            data
        )

        # Compute global matrices
        C, K, F = CKF_withpt(conM, cordM, H, k, c, htc, coord_sys, Tamb)

        B_n = (C - dt * (1 - theta) * K).dot(T_temp) + dt * (1 - theta) * F

        ss = 0
        epmax = 1E-6

        # --- Nonlinear iteration loop ---
        while ep > epmax:

            ss += 1

            H, k, c, _, _, _, _, _, _, _ = ThermalPropEnthalpy(
                T_temp,
                XM[:, ts - 1],
                XB[:, ts - 1],
                XP[:, ts - 1],
                XF[:, ts - 1],
                Sf, Sp, Sb,
                dt,
                data
            )

            _, K_temp, F_temp = CKF_withpt(conM, cordM, H, k, c, htc, coord_sys, Tamb)

            A = dt * theta * K_temp + C
            B = B_n + dt * theta * F_temp

            # Solve system
            T_new = np.linalg.pinv(A, rtol=1e-4).dot(B)

            T_new = np.maximum(T_new, Tamb)
            # Convert column vector to 1-D vector
            T_new = np.asarray(T_new).ravel()

            ep = np.sqrt(
                np.sum((T_new - T_temp.ravel()) ** 2)
                /
                np.sum(T_new ** 2)
            )
            # Ensure real temperature
            if not np.isrealobj(T_new):
                T_new = np.real(T_new)
            T_temp = T_new.copy()

            # loosen epmax if converging too slowly
            # Only relax tolerance if we've tried enough
            if ss > 10:
                epmax *= 10
                if epmax >= 1e-1:
                    break

        T = np.column_stack((T, T_temp))
        xa = in_aus - xm - xb - xp - xf
        XA = np.column_stack((XA, xa))
        XM = np.column_stack((XM, xm))
        XB = np.column_stack((XB, xb))
        XP = np.column_stack((XP, xp))
        XF = np.column_stack((XF, xf))

        print(f"\ntime step = {ts}")

        if T_temp[0] <= Tamb + 1:
            T_temp[:] = Tamb
            T[:, ts] = T_temp
            break

        if ts % 1000 == 0:
            print(
                f"Step={ts}, "
                f"Tcenter={T_temp[0]:.3f}, "
                f"Tsurface={T_temp[-1]:.3f}"
            )

    # --- Save temperature matrix T as TXT with tab delimiter ---
    S = np.column_stack((Sf, Sp, Sb))
    T_path = os.path.join(output_folder, 'temperature_matrix.txt').replace('\\', '/')
    xau_path = os.path.join(output_folder, 'Austenite_matrix.txt').replace('\\', '/')
    xfe_path = os.path.join(output_folder, 'Ferrite_matrix.txt').replace('\\', '/')
    xpe_path = os.path.join(output_folder, 'Pearlite_matrix.txt').replace('\\', '/')
    xba_path = os.path.join(output_folder, 'Bainite_matrix.txt').replace('\\', '/')
    xma_path = os.path.join(output_folder, 'Martensite_matrix.txt').replace('\\', '/')
    ss_path = os.path.join(output_folder, 'Schiel_Sum_matrix.txt').replace('\\', '/')
    np.savetxt(T_path, np.round(T, 2), delimiter='\t')
    np.savetxt(xau_path, np.round(XA, 4), delimiter='\t')
    np.savetxt(xfe_path, np.round(XF, 4), delimiter='\t')
    np.savetxt(xpe_path, np.round(XP, 4), delimiter='\t')
    np.savetxt(xba_path, np.round(XB, 4), delimiter='\t')
    np.savetxt(xma_path, np.round(XM, 4), delimiter='\t')
    np.savetxt(ss_path, np.round(S, 4), delimiter='\t')

    print(f"Results saved to {output_folder}")
