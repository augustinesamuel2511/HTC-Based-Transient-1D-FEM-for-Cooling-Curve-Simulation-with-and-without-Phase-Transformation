import numpy as np
import os


def femax_calc(c, cr, mn, mo, ni, si, output_folder):
    def Ce_from_composition(Cr, Mn, Mo, Ni, Si):
        """
        Solve Ae3(C) = Ae1(C) for C (Ce) given Cr, Mn, Mo, Ni, Si.
        Returns all real roots and the 'physical' Ce (0–1.5 wt.%).
        """
        # --- Ae3 coefficients: Ae3 = a2*C^2 + a1*C + a0 ---
        a2 = 90.91
        a1 = (-275.89
              + 16.45 * Cr
              - 29.96 * Mn
              - 10.80 * Mo)

        a0 = (
                883.49
                - 12.26 * Cr
                + 8.49 * Mo
                - 25.56 * Ni
                + 1.45 * Mn * Ni
                + 0.76 * Ni ** 2
                + 13.53 * Si
                - 3.47 * Mn * Si
        )

        # --- Ae1 coefficients: Ae1 = b1*C + b0 ---
        b1 = (-1.03 * Cr
              + 0.91 * Mn
              + 0.53 * Ni)

        b0 = (
                727.37 + 13.40 * Cr - 16.72 * Mn + 6.18 * Cr * Mn - 0.64 * Mn ** 2
                + 3.14 * Mo + 1.86 * Cr * Mo - 0.73 * Mn * Mo
                - 13.66 * Ni + 1.11 * Cr * Ni - 2.28 * Mn * Ni - 0.24 * Ni ** 2
                + 6.34 * Si - 8.88 * Cr * Si - 2.34 * Mn * Si + 11.98 * Si ** 2
        )

        # Equate Ae3(C) = Ae1(C):
        # a2*C^2 + a1*C + a0 = b1*C + b0
        # → a2*C^2 + (a1 - b1)*C + (a0 - b0) = 0
        A = a2
        B = a1 - b1
        Cc = a0 - b0

        roots = np.roots([A, B, Cc])
        real_roots = [r.real for r in roots if abs(r.imag) < 1e-8]

        # Pick a physically meaningful Ce (say 0–1.5 wt.%)
        valid = [r for r in real_roots if 0.0 <= r <= 1.5]
        Ce = valid[0] if valid else None

        return Ce, real_roots

    def Femax(C, Ce):
        """
        Femax = (Ce - C) / (Ce - 0.02)
        """
        if Ce is None:
            raise ValueError("Ce is None; cannot compute Femax.")
        if abs(Ce - 0.02) < 1e-12:
            raise ZeroDivisionError("Ce too close to 0.02.")
        return (Ce - C) / (Ce - 0.02)

    ce, root = Ce_from_composition(cr, mn, mo, ni, si)
    print("All roots of Ae3 = Ae1:", root)
    print("Chosen Ce =", ce)

    Fe_max = Femax(c, ce)
    print(Fe_max)

    os.makedirs(output_folder, exist_ok=True)
    femax_file = os.path.join(output_folder, "equilibrium_ferrite.txt").replace('\\', '/')
    print("Saving Fe_max to:", femax_file)

    with open(femax_file, "w") as f:
        f.write(f"{Fe_max:.6f}\n")

    print("File written.")
