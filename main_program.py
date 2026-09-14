import tkinter as tk
from tkinter import filedialog
import pandas as pd
import numpy as np
import os
from discretizing_meshing import discretizing_meshing
from equilibrium_ferrite_calc import femax_calc
from FEM1D_NoPT import FEM1D
from FEM1D_WithPT import FEM1D_withpt
from plotting import GraphPlotter


def browse_directory():
    folder_selected = filedialog.askdirectory()
    if folder_selected:  # if user picked a folder
        dir_entry.delete(0, tk.END)  # clear previous text
        dir_entry.insert(0, folder_selected)  # insert new path


root = tk.Tk()
root.title("1D Transient Heat Conduction Equation")
root.geometry("1200x750")
root.minsize(1100, 700)

# Columns
root.grid_columnconfigure(0, weight=1)
root.grid_columnconfigure(1, weight=1)
root.grid_columnconfigure(2, weight=2)

# Rows
root.grid_rowconfigure(0, weight=1)
root.grid_rowconfigure(1, weight=1)

frame_1 = tk.LabelFrame(root, text="Output directory", padx=10, pady=10)
frame_1.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)  # Place frame in grid at row 0, col 0

# Label
tk.Label(frame_1, text="Output Directory:", font=("Arial", 12), padx=5, pady=5).grid(column=0, row=0)

# Entry field
dir_entry = tk.Entry(frame_1, width=25, font=("Arial", 12))
dir_entry.grid(column=1, row=0, pady=10)

# Browse button
browse_btn = tk.Button(frame_1, text="Browse", font=("Arial", 12), command=browse_directory)
browse_btn.grid(column=2, row=0, pady=10)


def get_output_folder():
    output_folder = dir_entry.get().strip().replace("\\", "/")

    if not output_folder:
        print("Please select an output directory.")
        return None

    os.makedirs(output_folder, exist_ok=True)

    return output_folder


# ================= Left Frame: Geometry and Mesh details =================
frame_2 = tk.LabelFrame(root, text="Geometry and Mesh details", padx=10, pady=10)
frame_2.grid(row=0, column=1, sticky="nsew", padx=10, pady=10)

# Coordinate system selection (radio buttons)
coord_var = tk.StringVar(value="None")  # Default value = None
tk.Label(frame_2, text="Coordinate System").grid(row=0, column=0, sticky="w")
tk.Radiobutton(frame_2, text="Cartesian", variable=coord_var, value="Cartesian").grid(row=1, column=0,
                                                                                      sticky="w")
tk.Radiobutton(frame_2, text="Cylindrical", variable=coord_var, value="Cylindrical").grid(row=2, column=0,
                                                                                          sticky="w")

# X/R dimension entry
tk.Label(frame_2, text="X/R dimension (mm)").grid(row=3, column=0, sticky="w")
x_dim = tk.Entry(frame_2)
x_dim.grid(row=3, column=1)

# X/R division entry
tk.Label(frame_2, text="X/R division").grid(row=4, column=0, sticky="w")
x_div = tk.Entry(frame_2)
x_div.grid(row=4, column=1)


def confirm_fn():
    output_folder = get_output_folder()

    if output_folder is None:
        return
    coord_file = os.path.join(output_folder, "coordinate_system.txt").replace('\\', '/')
    print(coord_file)
    with open(coord_file, "w") as f:
        f.write(coord_var.get() + "\n")
    try:
        dimension = float(x_dim.get())
        divisions = int(x_div.get())
    except ValueError:
        print("Enter valid numeric values for dimension and division.")
        return

    if dimension <= 0:
        print("Dimension must be greater than zero.")
        return

    if divisions <= 0:
        print("Division must be greater than zero.")
        return

    discretizing_meshing(
        dimension,
        divisions,
        output_folder
    )


# Confirm button for Geometry section
tk.Button(frame_2, text="Confirm", command=confirm_fn).grid(row=8, column=0, columnspan=2, pady=10)

# ================= Material details =================

frame_3 = tk.LabelFrame(
    root,
    text="Material details",
    padx=10,
    pady=10
)

frame_3.grid(
    row=0,
    column=2,
    rowspan=2,
    sticky="nsew"
)

phase_var = tk.StringVar(value="None")

tk.Radiobutton(frame_3, text="Without phase transformation",
               variable=phase_var, value="noPT",
               command=lambda: show_material_frame()).grid(row=0, column=0, sticky="w")

tk.Radiobutton(frame_3, text="With phase transformation",
               variable=phase_var, value="withPT",
               command=lambda: show_material_frame()).grid(row=1, column=0, sticky="w")

# ======= WITHOUT PHASE TRANSFORMATION =======
frame_NPT = tk.Frame(frame_3, highlightbackground="gray", highlightthickness=1, padx=10, pady=10)


# Helper function to add each material property section
def add_material_section(parent, title, start_row):
    # Title
    tk.Label(parent, text=title).grid(row=start_row, column=0, sticky="w")

    entries = {}  # dict to hold entry widgets

    # Create labels and entry boxes
    for i, lbl in enumerate(["a1", "a2", "a3", "Constant"]):
        tk.Label(parent, text=lbl).grid(row=start_row + i + 1, column=0, sticky="w")
        ent = tk.Entry(parent)
        ent.grid(row=start_row + i + 1, column=1)
        entries[lbl] = ent  # store entry by label name

    return start_row + 5, entries


# Add sections for Density, Thermal Conductivity, Specific Heat
row = 0
row, density_entries = add_material_section(frame_NPT, "Density (Kg/m^3)", row)
row, k_entries = add_material_section(frame_NPT, "Thermal Conductivity (W/m K)", row)
row, cp_entries = add_material_section(frame_NPT, "Specific Heat (J/Kg K)", row)

# ===================== WITH PHASE TRANSFORMATION FRAME =====================

frame_PT = tk.Frame(
    frame_3,
    highlightbackground="gray",
    highlightthickness=1,
    pady=10
)

tk.Label(frame_PT, text="Temperature and Phase Dependent Thermal Properties",
         font=("Arial", 10, "bold")).grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 5))

tk.Label(frame_PT, text="Temperature in degree Celsius",
         font=("Arial", 10, "italic")).grid(row=1, column=0, columnspan=3, sticky="w")


# ------- Helper function for browse & entries -------
def browse_excel(target_entry):
    filepath = filedialog.askopenfilename(
        title="Select Excel/CSV/TXT file",
        filetypes=[("Excel/CSV/TXT", "*.xlsx *.xls *.csv *.txt"), ("All files", "*.*")]
    )
    if filepath:
        target_entry.delete(0, tk.END)
        target_entry.insert(0, filepath)


# Row counter for PT properties
r = 2

# ------- Thermal Conductivity -------
tk.Label(frame_PT, text="Thermal conductivity (W/mK)").grid(row=r, column=0, sticky="e")
tc_PT_entry = tk.Entry(frame_PT, width=25)
tc_PT_entry.grid(row=r, column=1)
tk.Button(frame_PT, text="Browse Excel file",
          command=lambda: browse_excel(tc_PT_entry)).grid(row=r, column=2, padx=2)
r += 1

# ------- Density -------
tk.Label(frame_PT, text="Density (Kg/m3)").grid(row=r, column=0, sticky="e")
rho_PT_entry = tk.Entry(frame_PT, width=25)
rho_PT_entry.grid(row=r, column=1)
tk.Button(frame_PT, text="Browse Excel file",
          command=lambda: browse_excel(rho_PT_entry)).grid(row=r, column=2, padx=2)
r += 1

# ------- Specific Heat -------
tk.Label(frame_PT, text="Specific heat (J/KgK)").grid(row=r, column=0, sticky="e")
cp_PT_entry = tk.Entry(frame_PT, width=25)
cp_PT_entry.grid(row=r, column=1)
tk.Button(frame_PT, text="Browse Excel file",
          command=lambda: browse_excel(cp_PT_entry)).grid(row=r, column=2, padx=2)
r += 1

# ------- Latent Heat -------
tk.Label(frame_PT, text="Latent heat (J/Kg)").grid(row=r, column=0, sticky="e")
Lh_PT_entry = tk.Entry(frame_PT, width=25)
Lh_PT_entry.grid(row=r, column=1)
tk.Button(frame_PT, text="Browse Excel file",
          command=lambda: browse_excel(Lh_PT_entry)).grid(row=r, column=2, padx=2)
r += 1

# ------- TTT Diagram -------
tk.Label(frame_PT, text="Load TTT Diagram").grid(row=r, column=0, sticky="e")
ttt_PT_entry = tk.Entry(frame_PT, width=25)
ttt_PT_entry.grid(row=r, column=1)
tk.Button(frame_PT, text="Browse Excel file",
          command=lambda: browse_excel(ttt_PT_entry)).grid(row=r, column=2, padx=2)
r += 1

# -------- Critical Temperatures --------
tk.Label(frame_PT, text="Critical Temperatures in degree Celsius(Put Ae3=0 if no ferrite transformation)",
         font=("Arial", 10, "bold")).grid(row=r, column=0, columnspan=2, sticky="w", pady=(10, 5))
r += 1

tk.Label(frame_PT, text="Ae3").grid(row=r, column=0, sticky="e")
Ae3_entry = tk.Entry(frame_PT, width=15)
Ae3_entry.grid(row=r, column=1, sticky="w")
r += 1

tk.Label(frame_PT, text="Ae1").grid(row=r, column=0, sticky="e")
Ae1_entry = tk.Entry(frame_PT, width=15)
Ae1_entry.grid(row=r, column=1, sticky="w")
r += 1

tk.Label(frame_PT, text="Bs").grid(row=r, column=0, sticky="e")
Bs_entry = tk.Entry(frame_PT, width=15)
Bs_entry.grid(row=r, column=1, sticky="w")
r += 1

tk.Label(frame_PT, text="Ms").grid(row=r, column=0, sticky="e")
Ms_entry = tk.Entry(frame_PT, width=15)
Ms_entry.grid(row=r, column=1, sticky="w")
r += 1

tk.Label(frame_PT, text="Initial Austenite Phase Fraction (0-1)").grid(row=r, column=0, sticky="e")
in_aus_entry = tk.Entry(frame_PT, width=15)
in_aus_entry.grid(row=r, column=1, sticky="w")
r += 1

tk.Label(frame_PT,
         text="Steel Composition in Wt% for equilibrium Ferrite Calculation (Put 0.0 for non existing element)",
         font=("Arial", 10, "bold")).grid(row=r, column=0, columnspan=2, sticky="w", pady=(10, 5))
tk.Label(frame_PT, text="C").grid(row=r + 1, column=0, sticky="e")
c_entry = tk.Entry(frame_PT, width=15)
c_entry.grid(row=r + 1, column=1, sticky="w")
tk.Label(frame_PT, text="Cr").grid(row=r + 2, column=0, sticky="e")
cr_entry = tk.Entry(frame_PT, width=15)
cr_entry.grid(row=r + 2, column=1, sticky="w")
tk.Label(frame_PT, text="Mn").grid(row=r + 3, column=0, sticky="e")
mn_entry = tk.Entry(frame_PT, width=15)
mn_entry.grid(row=r + 3, column=1, sticky="w")
tk.Label(frame_PT, text="Mo").grid(row=r + 4, column=0, sticky="e")
mo_entry = tk.Entry(frame_PT, width=15)
mo_entry.grid(row=r + 4, column=1, sticky="w")
tk.Label(frame_PT, text="Ni").grid(row=r + 5, column=0, sticky="e")
ni_entry = tk.Entry(frame_PT, width=15)
ni_entry.grid(row=r + 5, column=1, sticky="w")
tk.Label(frame_PT, text="Si").grid(row=r + 6, column=0, sticky="e")
si_entry = tk.Entry(frame_PT, width=15)
si_entry.grid(row=r + 6, column=1, sticky="w")


def show_material_frame():
    frame_NPT.grid_forget()
    frame_PT.grid_forget()

    if phase_var.get() == "noPT":
        frame_NPT.grid(row=3, column=0, sticky="nw", padx=5, pady=5)
        microstructure_btn.grid_remove()

    elif phase_var.get() == "withPT":
        frame_PT.grid(row=3, column=0, sticky="nw", padx=5, pady=5)
        microstructure_btn.grid(row=3, column=0, pady=20)


def load_property_file(filepath):
    """Reads an Excel/CSV/TXT file and returns numeric data as a 2D list."""

    filepath = filepath.strip()
    if not filepath:
        return None

    try:
        ext = filepath.lower()

        if ext.endswith((".xlsx", ".xls")):
            df = pd.read_excel(filepath, header=None)

        elif ext.endswith(".csv"):
            df = pd.read_csv(filepath, header=None)

        elif ext.endswith(".txt"):
            df = pd.read_csv(filepath, header=None, sep=None, engine="python")

        else:
            raise ValueError("Unsupported file type for: " + filepath)

        # Convert to numeric only
        df = df.apply(pd.to_numeric, errors="coerce").dropna(axis=1, how="all").dropna(axis=0, how="all")

        return df.values.tolist()  # return as list of lists

    except Exception as e:
        print("Error reading property file:", filepath, e)
        return None


def save_material_properties():
    output_folder = get_output_folder()

    if output_folder is None:
        return
    os.makedirs(output_folder, exist_ok=True)
    pt_file = os.path.join(output_folder, "PT_file.txt").replace('\\', '/')
    with open(pt_file, "w") as f:
        f.write(phase_var.get())

    if phase_var.get() == "noPT":
        mat_file = os.path.join(output_folder, "material_properties.txt").replace('\\', '/')

        # Build 3x4 matrix: rows = properties, cols = [a1, a2, a3, Constant]
        mat = np.zeros((3, 4), dtype=float)

        # Fill from entries
        mat[0, :] = [float(density_entries[k].get() or 0) for k in ["a1", "a2", "a3", "Constant"]]
        mat[1, :] = [float(k_entries[k].get() or 0) for k in ["a1", "a2", "a3", "Constant"]]
        mat[2, :] = [float(cp_entries[k].get() or 0) for k in ["a1", "a2", "a3", "Constant"]]

        # Save as text file
        np.savetxt(mat_file, mat, fmt="%.6f")
        print("Material properties saved as 3x4 matrix in:", mat_file)

    elif phase_var.get() == "withPT":
        # --- 1. Read ALL property Excel/CSV/TXT files ---
        tc_data = load_property_file(tc_PT_entry.get())
        rho_data = load_property_file(rho_PT_entry.get())
        cp_data = load_property_file(cp_PT_entry.get())
        lh_data = load_property_file(Lh_PT_entry.get())
        ttt_data = load_property_file(ttt_PT_entry.get())

        # --- 2. Save each into its own TXT file ---
        files_and_data = {
            "thermal_conductivity_PT.txt": tc_data,
            "density_PT.txt": rho_data,
            "specific_heat_PT.txt": cp_data,
            "latent_heat_PT.txt": lh_data,
            "ttt_diagram.txt": ttt_data,
        }

        for filename, data in files_and_data.items():
            out_file = os.path.join(output_folder, filename).replace('\\', '/')

            if data is None:
                with open(out_file, "w") as f:
                    f.write("ERROR: Could not read data\n")
            else:
                np.savetxt(out_file, np.array(data), fmt="%.6f")

            print("Saved:", out_file)

        # --- Critical Temperatures (ROW VECTOR) ---
        crit_file = os.path.join(output_folder, "critical_temperatures.txt").replace('\\', '/')
        in_aus_file = os.path.join(output_folder, "initial_austenite.txt").replace('\\', '/')
        Ae3_val = Ae3_entry.get()
        Ae1_val = Ae1_entry.get()
        Bs_val = Bs_entry.get()
        Ms_val = Ms_entry.get()
        in_aus_val = in_aus_entry.get()

        with open(crit_file, "w") as f:
            f.write(f"{Ae3_val} {Ae1_val} {Bs_val} {Ms_val}\n")

        with open(in_aus_file, "w") as f:
            f.write(in_aus_val)

        c_val = float(c_entry.get())
        cr_val = float(cr_entry.get())
        mn_val = float(mn_entry.get())
        mo_val = float(mo_entry.get())
        ni_val = float(ni_entry.get())
        si_val = float(si_entry.get())
        femax_calc(c=c_val, cr=cr_val, mn=mn_val, mo=mo_val, ni=ni_val, si=si_val, output_folder=output_folder)
        print("Saved WITH-PT material files:")


# Confirm button for Material section
tk.Button(frame_3, text="Confirm", command=save_material_properties).grid(row=row, column=0, columnspan=2,

                                                                          pady=10)
frame_4 = tk.LabelFrame(root, text="Solver details", padx=10, pady=10)
frame_4.grid(row=1, column=0, sticky="nsew")
# HTC input
tk.Label(frame_4, text="Boundary Heat Transfer Coefficient (W/m^2K)").grid(row=0, column=0, sticky="w")
htc = tk.Entry(frame_4)
htc.grid(row=0, column=1)

# Initial Temperature
tk.Label(frame_4, text="Initial Temperature (degC)").grid(row=1, column=0, sticky="w")
Tin = tk.Entry(frame_4)
Tin.grid(row=1, column=1)

# Ambient temperature input
tk.Label(frame_4, text="Ambient Temperature (degC)").grid(row=2, column=0, sticky="w")
Tamb = tk.Entry(frame_4)
Tamb.grid(row=2, column=1)

# time interval dt
tk.Label(frame_4, text="time interval, dt (s)").grid(row=3, column=0, sticky="w")
dt = tk.Entry(frame_4)
dt.grid(row=3, column=1)


def solve():
    output_folder = get_output_folder()

    if output_folder is None:
        return

    try:
        Tin_value = float(Tin.get())
        Tamb_value = float(Tamb.get())
        htc_value = float(htc.get())
        dt_value = float(dt.get())
    except ValueError:
        print("Please enter valid solver values.")
        return
    dt_path = os.path.join(output_folder, 'dt.txt').replace('\\', '/')
    np.savetxt(dt_path, np.array([np.round(dt_value, 2)])) 
    pt_file = os.path.join(
        output_folder,
        "PT_file.txt"
    )

    if not os.path.exists(pt_file):
        print("Please confirm material properties first.")
        return

    with open(pt_file, "r") as f:
        pt_flag = f.read().strip()

    if pt_flag == "noPT":

        FEM1D(
            Tin_value,
            Tamb_value,
            htc_value,
            dt_value,
            output_folder
        )

    elif pt_flag == "withPT":

        FEM1D_withpt(
            Tin_value,
            Tamb_value,
            htc_value,
            dt_value,
            output_folder
        )

    else:
        print("Error in material property selection.")


# Solve button
tk.Button(frame_4, text="Solve", command=solve).grid(row=4, column=0, columnspan=2, pady=20)

frame_5 = tk.LabelFrame(root, text="Post Process")
frame_5.grid(row=1, column=1, sticky="nsew", padx=10, pady=10)
# Instruction label
tk.Label(frame_5, text="Provide X/R locations separated by comma").grid(row=0, column=0, columnspan=2,
                                                                        sticky="w", pady=5)

# Entry for X/R locations
tk.Label(frame_5, text="X/R location (mm)").grid(row=1, column=0, sticky="w")
node_loc = tk.Entry(frame_5, width=40)
node_loc.grid(row=1, column=1, padx=5)


def temp_plot():
    output_folder = get_output_folder()

    if output_folder is None:
        return
    DT_file = os.path.join(output_folder, "dt.txt").replace('\\', '/')
    DT = np.loadtxt(DT_file)
    loc = node_loc.get()
    loc_lst = [float(x.strip()) for x in loc.split(",")]
    os.makedirs(output_folder, exist_ok=True)

    cord_file = os.path.join(output_folder, "cordM.txt").replace('\\', '/')
    cordM = np.loadtxt(cord_file)

    temp_file = os.path.join(output_folder, "temperature_matrix.txt").replace('\\', '/')
    temp_mat = np.loadtxt(temp_file)  # shape: (n_coord, n_time)

    time = np.arange(temp_mat.shape[1]) * DT
    # Find indices for each requested location
    index = {}
    for x in loc_lst:
        exact = np.where(np.isclose(cordM, x))[0]
        if exact.size > 0:
            index[x] = int(exact[0])
        elif x < cordM[0] or x > cordM[-1]:
            raise ValueError(f"{x} is outside the coordinate range.")
        else:
            i = np.searchsorted(cordM, x)
            index[x] = [i - 1, i]

    # Extract or interpolate temperatures
    temp_curves = []
    for x, idx in index.items():
        if isinstance(idx, int):  # exact match
            temp = temp_mat[idx, :]  # row → curve across time
        else:  # interpolation between two coords
            i1, i2 = idx
            x1, x2 = cordM[i1], cordM[i2]
            frac = (x - x1) / (x2 - x1)
            temp = temp_mat[i1, :] + frac * (temp_mat[i2, :] - temp_mat[i1, :])
        temp_curves.append(temp)

    # Build matrix: time + all curves (transpose temp_curves to match time length)
    temp_curve_matrix = np.column_stack([time] + temp_curves)

    # Plot
    gp1 = GraphPlotter("Temperature Curves", "Time (s)", "Temperature (°C)")
    gp1.plot_single_multi_y(x=time, y_list=temp_curves, labels=[str(x) for x in index.keys()])
    gp1.show()

    # Save to CSV
    extract = os.path.join(output_folder, "Temperature Curve Extract.csv").replace('\\', '/')
    header = "time," + ",".join([f"T(x={x})" for x in index.keys()])
    np.savetxt(extract, temp_curve_matrix, delimiter=",", header=header, comments="")


# Buttons: Plot and Extract
tk.Button(frame_5, text="Plot & Extract", command=temp_plot).grid(row=2, column=0, pady=20)


def microstructure_plot():
    output_folder = get_output_folder()

    if output_folder is None:
        return
    os.makedirs(output_folder, exist_ok=True)

    cord_file = os.path.join(output_folder, "cordM.txt").replace('\\', '/')
    xa_file = os.path.join(output_folder, "Austenite_matrix.txt").replace('\\', '/')
    xf_file = os.path.join(output_folder, "Ferrite_matrix.txt").replace('\\', '/')
    xp_file = os.path.join(output_folder, "Pearlite_matrix.txt").replace('\\', '/')
    xb_file = os.path.join(output_folder, "Bainite_matrix.txt").replace('\\', '/')
    xm_file = os.path.join(output_folder, "Martensite_matrix.txt").replace('\\', '/')
    cordM = np.loadtxt(cord_file)
    xa = np.loadtxt(xa_file)[:, -1]
    xf = np.loadtxt(xf_file)[:, -1]
    xp = np.loadtxt(xp_file)[:, -1]
    xb = np.loadtxt(xb_file)[:, -1]
    xm = np.loadtxt(xm_file)[:, -1]
    y_file = np.column_stack([xa, xf, xp, xb, xm])
    gp2 = GraphPlotter("Phase Distribution", "Distance (mm)", "Phase Fraction")
    gp2.plot_single_multi_y(x=cordM, y_list=y_file.T,
                            labels=['Austenite', 'Ferrite', 'Pearlite', 'Bainite', 'Martensite'])
    gp2.show()


microstructure_btn = tk.Button(frame_5, text="Plot Microstructure Distribution", command=microstructure_plot)

microstructure_btn.grid(row=3, column=0, pady=20)

microstructure_btn.grid_remove()

# Back button (bottom right corner)
tk.Button(frame_5, text="Close", command=root.destroy).grid(row=4, column=0, sticky="e", pady=10)

root.mainloop()
