import fem_toolbox.femsolver as ft
from fem_toolbox.elements import rotation_matrix_bar_2d
import numpy as np
import matplotlib.pyplot as plt
import inspect

def eval_stress(k_local_func, R_func, u, fem_elements, fem_nodes, element_sections, material, section_shape, ndof_per_node):
    """
    Computes stress in each element.

    Args:
        u: Global displacement vector.
        fem_elements: (n_elements, 2) connectivity matrix.
        fem_nodes: (n_nodes, 2) node coordinates.
        element_sections: (n_elements, 2) array of [A, I] for each element.
        section_shape: str - "rectangle" or "circular"
        material: {"E" = E, "rho" = rho}
        k_local_func: function returning 6x6 stiffness matrix in local system (so theoretically you could use a different function than k_beam2d)
        R_func: function returning rotation matrix.

    Returns:
        stress_max: max absolute stress (tension or compression) per element
        stress_axial: N / A
        stress_bending: M * c / I
        stress_shear: estimated max shear stress (section dependent)
        vonMises_stress
    """
    E = material['E']
    stress_max = []
    stress_axial = []
    stress_bending = []
    stress_shear = []
    von_mises_stress = []
    Nvec = []
    Vvec = []
    Mvec = []

    # Inspect the signature of the provided local stiffness function
    sig = inspect.signature(k_local_func)
    param_names = list(sig.parameters.keys())

    for e, (n1, n2) in enumerate(fem_elements):
        
        # get node coordinates for each element
        x1, y1 = fem_nodes[n1]
        x2, y2 = fem_nodes[n2]

        # get geometry params for each element
        A, I = element_sections[e]['A'], element_sections[e]['I']
        L = np.hypot(x2 - x1, y2 - y1)

        # get height based on the cross section shape
        match section_shape:
            case "rectangle":
                h = np.sqrt(12 * I / A)
                c = h / 2
            case "circular":
                c = np.sqrt(A / np.pi)
            case _ :
                return "Section shape not supported. Only possibilities are 'rectangle' and 'circular'"

        # DOF mapping
        flag = False
        if ndof_per_node == 1 and R_func is not None:
            ndof_per_node = 2 # correct for 2d truss structures
            flag = True

        dof_map = ft.dof_map(ndof_per_node, n1, n2)
        # displacement of the 2 nodes of the element in the global system
        u_e_global = u[dof_map]


        # Prepare all possible arguments
        available_args = {
            'E': E,
            'A': A,
            'I': I,
            'L': L
        }


        # Filter only the arguments needed for k_local_func
        args_to_pass = [available_args[name] for name in param_names if name in available_args]

        # rotation matrix from local to global
        if flag:
            ndof_per_node = 1

        if R_func == None:
            if ndof_per_node == 1:
                R = np.eye(2)
            else:
                R = np.eye(4)
        else:
            if ndof_per_node == 2:  # need to cut away raws and columns 3 and 6
                R = np.delete(R_func(x1,y1,x2,y2), [2, 5], axis=0)  # cut rows
                R = np.delete(R, [2, 5], axis=1)
            
            elif ndof_per_node == 1:
                R = rotation_matrix_bar_2d(x1, y1, x2, y2)     # 4x4

            else:
                R = R_func(x1, y1, x2, y2)

        # rotate to local coordinates
        u_e_local = R @ u_e_global

        # Call with only the required arguments
        k_local = k_local_func(*args_to_pass)
        f_local = k_local @ u_e_local

        #print(f"u local: {u_e_local}\nf_local: {f_local}")

        match section_shape:
            case "rectangle":
                coeff = 3/2
            case "circular":
                coeff = 4/3
            case _ :
                return "Section shape not supported. Only possibilities are 'rectangle' and 'circular'"
            
        if flag:
            ndof_per_node = 1

        if ndof_per_node == 3:
            # axial force
            delta_L = u_e_local[3] - u_e_local[0]  # u2_local - u1_local (in axial dir)
            strain = delta_L / L
            N = E * A * strain

            # shear force
            V = f_local[1]
            # bending moment 
            M = f_local[2]

            # Stress components
            sigma_axial = N / A
            sigma_bend = - M * c / I    # negative because positive M compresses the top

            tau_xy = coeff * V / A

            # Total max tension/compression stress (superposition)

            sigma_top = sigma_axial + sigma_bend 
            sigma_bottom = sigma_axial - sigma_bend 

            # take the max, preserving the sign
            sigma_max = sigma_top if abs(sigma_top) >= abs(sigma_bottom) else sigma_bottom


            # vonMises stress
            sigma_von_mises = np.sqrt(sigma_max**2 + 3*tau_xy**2)
        
            stress_max.append(sigma_max)
            stress_axial.append(sigma_axial)
            stress_bending.append(sigma_bend)
            stress_shear.append(tau_xy)
            von_mises_stress.append(sigma_von_mises)

            Nvec.append(N)
            Vvec.append(V)
            Mvec.append(M)

            internal_actions = np.column_stack((Nvec, Vvec, Mvec))

        if ndof_per_node == 2:
            # shear force
            V = f_local[0]
            # bending moment 
            M = f_local[1]

            # Stress components
            sigma_bend =  M * c / I  

            tau_xy = coeff * V / A

            # Total max tension/compression stress (superposition)
            sigma_max = sigma_bend 

            # vonMises stress
            sigma_von_mises = np.sqrt(sigma_max**2 + 3*tau_xy**2)

            stress_max.append(sigma_max)
            stress_bending.append(sigma_bend)
            stress_shear.append(tau_xy)
            von_mises_stress.append(sigma_von_mises)

            Vvec.append(V)
            Mvec.append(M)

            internal_actions = np.column_stack((Vvec, Mvec))

        if ndof_per_node == 1:
            # axial force
            delta_L = u_e_local[1] - u_e_local[0]  # u2_local - u1_local (in axial dir)
            strain = delta_L / L
            N = E * A * strain
#            N = f_local[0]

            # Stress components
            sigma_axial = N / A
            # Total max tension/compression stress (superposition)
            sigma_max = sigma_axial 

            # vonMises stress
            sigma_von_mises = sigma_axial 


            stress_max.append(sigma_max)
            stress_axial.append(sigma_axial)
            von_mises_stress.append(sigma_von_mises)

            Nvec.append(N)

            internal_actions = np.column_stack((Nvec))


    return (
        np.array(stress_max), 
        np.array(stress_axial), 
        np.array(stress_bending),
        np.array(stress_shear),
        np.array(von_mises_stress),
        internal_actions
    )



def compute_reaction_forces(K_global, U, bc_nodes, bc_dofs, ndof_per_node):
    """
    Compute the reaction forces at constrained DOFs.

    Parameters:
        K_global: Full global stiffness matrix (N_dof x N_dof)
        U: Full global displacement vector (including constrained DOFs)
        bc_nodes: List of nodes with prescribed BCs
        bc_dofs:  List of DOFs per node (e.g. 0 for ux, 1 for uy, 2 for theta)

    Returns:
        reactions: Array of reaction forces at prescribed DOFs
        constrained_dof_ids: List of constrained DOF indices used
    """
    flag = False
    # understand if 2d truss structure
    if ndof_per_node == 1 and any(i >= 1 for i in bc_dofs):
        ndof_per_node = 2
        flag = True

    # Map (node, dof) to global DOF indices
    constrained_dof_ids = [node * ndof_per_node + dof for node, dof in zip(bc_nodes, bc_dofs)]

    # corrections if using beam or bar elements
    if ndof_per_node == 2 and not flag:
        # shift indexes
        bc_dofs = np.array(bc_dofs)
        bc_dofs -= 1
        constrained_dof_ids = [node * ndof_per_node + dof for node, dof in zip(bc_nodes, bc_dofs)]


    # Compute internal force vector from stiffness * displacement
    F_internal = K_global @ U

    # Reactions are internal forces at constrained DOFs
    reactions = F_internal[constrained_dof_ids]
    
    for dof, r in zip(constrained_dof_ids, reactions):
        print(f"DOF {dof}: Reaction force = {r:.2f} N")

    return reactions, constrained_dof_ids


def plot_internal_actions_2D(
    beam_connectivity,
    fem_nodes,
    fem_elements,
    internal_actions,
    elements_per_beam=5
):
    """
    Plot internal action diagrams (Axial, Shear, Moment) per beam.
    
    Each beam is shown with:
    - Black dashed line as the physical beam
    - Internal force/moment diagram plotted on top

    Parameters:
        node_coords_beam: (n_beam_nodes, 2) array of node coordinates (coarse geometry).
        beam_connectivity: (n_beams, 2) array of beam connectivity (coarse definition).
        fem_nodes: (n_fem_nodes, 2) array of FEM node coordinates.
        fem_elements: (n_elements, 2) array of FEM element connectivity.
        internal_actions: (n_elements, 3) array of [N, V, M] values per element.
        elements_per_beam: number of FEM elements per beam (default = 5).
    """
    num_beams = beam_connectivity.shape[0]

    for beam_id in range(num_beams):
        el_start = beam_id * elements_per_beam
        el_end = el_start + elements_per_beam
        beam_elements = fem_elements[el_start:el_end]
        beam_actions = internal_actions[el_start:el_end]  # shape: (elements_per_beam, 3)

        arc_lengths = []
        N_vals, V_vals, M_vals = [], [], []

        base_node = fem_nodes[beam_elements[0][0]]

        for i in range(elements_per_beam):
            n1_idx = beam_elements[i][0]
            pt = fem_nodes[n1_idx]
            arc_len = np.linalg.norm(pt - base_node)

            arc_lengths.append(arc_len)
            
            try:
                N_vals.append(beam_actions[i][0])
                V_vals.append(beam_actions[i][1])
                M_vals.append(beam_actions[i][2])
            
            except IndexError as e:
                raise IndexError(
                    "Indexing error: this may occur if you are trying to use this function for bar or beam elements. If you are studying a truss structure, the only internal actions are axial. If you are trying to study a structure made by simple beam elements (not 2D beams), you may instead want to use 2Dbeam elements (3 dofs)   "
                )

        # Add last node position and extrapolate values
        pt_last = fem_nodes[beam_elements[-1][1]]
        arc_last = np.linalg.norm(pt_last - base_node)
        arc_lengths.append(arc_last)

        def extrapolate(values):
            return values[-1] + (values[-1] - values[-2])

        N_vals.append(extrapolate(N_vals))
        V_vals.append(extrapolate(V_vals))
        M_vals.append(extrapolate(M_vals))

        # Plotting
        fig, axs = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
        fig.suptitle(f'Internal Actions – Beam {beam_id + 1}', fontsize=14)

        margin = 0.05 * (arc_lengths[-1] - arc_lengths[0])
        xlim = [arc_lengths[0] - margin, arc_lengths[-1] + margin]

        # Beam baseline for reference
        for ax in axs:
            ax.plot([arc_lengths[0], arc_lengths[-1]], [0, 0], 'k--', linewidth=1.5, label='Beam centerline')
            ax.set_xlim(xlim)
            ax.legend(loc="upper right")
            ax.grid(True)

        axs[0].plot(arc_lengths, N_vals, '-o', color='tab:blue', linewidth=2.5, label='Axial Force')
        axs[0].set_ylabel('Axial [N]')
        axs[0].legend()

        axs[1].plot(arc_lengths, V_vals, '-o', color='tab:orange', linewidth=2.5, label='Shear Force')
        axs[1].set_ylabel('Shear [N]')
        axs[1].legend()

        axs[2].plot(arc_lengths, M_vals, '-o', color='tab:green', linewidth=2.5, label='Bending Moment')
        axs[2].set_ylabel('Moment [Nm]')
        axs[2].set_xlabel('Length along beam [m]')
        axs[2].legend()

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.show()


def plot_internal_actions_1D(fem_nodes, fem_elements, internal_actions, dof_per_node):
    """
    Plots internal actions (axial, shear, moment) for 1D horizontal geometries.

    Parameters:
        fem_nodes : (n_nodes, 2) array of node coordinates (only x is used)
        fem_elements : (n_elements, 2) array of node indices
        internal_actions : (n_elements, N) array of internal actions per element
            - N = 1 for axial (1 DOF)
            - N = 2 for shear + moment (2 DOFs)
        dof_per_node : int, either 1 or 2
    """
    x_coords = fem_nodes[:, 0]  # assume 1D geometry along X

    x_vals = []
    axial_vals = []
    shear_vals = []
    moment_vals = []

    for i, (n1, n2) in enumerate(fem_elements):
        x1 = x_coords[n1]
        x_vals.append(x1)

        if dof_per_node == 1:
            N = internal_actions[i][0]
            axial_vals.append(N)
        #if dof_per_node == 1:
        #    N = internal_actions[i] if np.isscalar(internal_actions[i]) else internal_actions[i][0]
        #    axial_vals.append(N)

        elif dof_per_node == 2:
            V, M = internal_actions[i]
            shear_vals.append(V)
            moment_vals.append(M)
        else:
            raise ValueError("Only dof_per_node = 1 or 2 supported in this 1D plotting function.")

    # Add final node (interpolated value)
    x_last = x_coords[fem_elements[-1][1]]
    x_vals.append(x_last)

    if dof_per_node == 1:
        # linear extrapolation for last point
        delta = axial_vals[-1] - axial_vals[-2]
        axial_vals.append(axial_vals[-1] + delta)
    else:
        delta_v = shear_vals[-1] - shear_vals[-2]
        delta_m = moment_vals[-1] - moment_vals[-2]
        shear_vals.append(shear_vals[-1] + delta_v)
        moment_vals.append(moment_vals[-1] + delta_m)

    # Create plot
    fig, axs = plt.subplots(2 if dof_per_node == 2 else 1, 1, figsize=(10, 4), sharex=True)
    if dof_per_node == 1:
        axs = [axs]  # make iterable

    if dof_per_node == 1:
        axs[0].plot(x_vals, axial_vals, 'r-', label='Axial force', linewidth=2)
        axs[0].plot(x_vals, [0]*len(x_vals), 'k--', linewidth=1.5)
        axs[0].plot(x_vals, axial_vals, 'ko', markersize=4)
        axs[0].set_ylabel("Axial N")
        axs[0].legend()
        axs[0].grid(True)

    else:
        axs[0].plot(x_vals, shear_vals, 'b-', label='Shear force', linewidth=2)
        axs[0].plot(x_vals, [0]*len(x_vals), 'k--', linewidth=1.5)
        axs[0].plot(x_vals, shear_vals, 'ko', markersize=4)
        axs[0].set_ylabel("Shear V")
        axs[0].legend()
        axs[0].grid(True)

        axs[1].plot(x_vals, moment_vals, 'g-', label='Bending moment', linewidth=2)
        axs[1].plot(x_vals, [0]*len(x_vals), 'k--', linewidth=1.5)
        axs[1].plot(x_vals, moment_vals, 'ko', markersize=4)
        axs[1].set_ylabel("Moment M")
        axs[1].set_xlabel("Position along beam (X)")
        axs[1].legend()
        axs[1].grid(True)

    plt.tight_layout()
    plt.show()


from matplotlib.collections import LineCollection
from matplotlib import colormaps

def plot_2D_loaded_structure(
    fem_nodes, fem_elements, u, stress_values, stress_type, ndof_per_node, scale=1.0, title=None, show_labels=False
):
    """
    Plots the undeformed and deformed structure with stress color map.

    Parameters:
        fem_nodes: (N, 2) array of node coordinates
        fem_elements: (E, 2) array of element connectivity
        u: global displacement vector 
        stress_values: (E,) array of stress values per element
        stress_type: str, one of 'von Mises', 'axial', 'bending', 'shear
        scale: float, deformation magnification factor
        title: plot title
        show_labels: if True, adds stress value labels at element midpoints
    """
    num_nodes = fem_nodes.shape[0]
    u_nodes = np.zeros_like(fem_nodes)

    # Extract displacements (ux, uy)
    for i in range(num_nodes):
        u_nodes[i, 0] = u[i * ndof_per_node]
        u_nodes[i, 1] = u[i * ndof_per_node + 1]

    # Apply scaling
    deformed_nodes = fem_nodes + scale * u_nodes

    # Prepare segments and colors
    lines_undeformed = []
    lines_deformed = []
    stress_colors = []

    for e, (n1, n2) in enumerate(fem_elements):
        p1, p2 = fem_nodes[n1], fem_nodes[n2]
        q1, q2 = deformed_nodes[n1], deformed_nodes[n2]

        lines_undeformed.append([p1, p2])
        lines_deformed.append([q1, q2])
        #stress_colors.append(von_mises_stress[e])
        stress_colors.append(stress_values[e])

    # Convert stress values to colormap
    stress_colors = np.array(stress_colors)
    cmap = colormaps['viridis']
    color_norm = (stress_colors - stress_colors.min()) / (stress_colors.max() - stress_colors.min() + 1e-9)
    color_mapped = cmap(color_norm)

    # Plotting
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.set_aspect('equal')

    # Undeformed (gray)
    undeformed = LineCollection(lines_undeformed, colors='gray', linestyle="--", linewidths=2, label='Undeformed')
    ax.add_collection(undeformed)

    # Deformed with stress color
    deformed = LineCollection(lines_deformed, colors=color_mapped, linewidths=3, label='Deformed')
    ax.add_collection(deformed)

    # Optional stress labels at element midpoints
    if show_labels:
        for e, (n1, n2) in enumerate(fem_elements):
            midpoint = (deformed_nodes[n1] + deformed_nodes[n2]) / 2
            ax.text(
                midpoint[0], midpoint[1],
                f"{stress_colors[e]:.1f}",
                color='black', fontsize=8, ha='center', va='center',
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2')
            )

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap)
    sm.set_array(stress_colors)
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label(f"{stress_type.capitalize()} Stress [MPa]")

    # Labels and axes
    if title is None:
        title = f"Deformation and {stress_type.capitalize()} Stress"
    ax.set_title(title)
    ax.set_xlabel("X [mm]")
    ax.set_ylabel("Y [mm]")
    ax.grid(True)
    ax.legend()

    # Set axis limits
    all_points = np.vstack((fem_nodes, deformed_nodes))
    x_min, x_max = all_points[:, 0].min(), all_points[:, 0].max()
    y_min, y_max = all_points[:, 1].min(), all_points[:, 1].max()
    pad_x = 0.1 * (x_max - x_min) if x_max != x_min else 1.0
    pad_y = 0.1 * (y_max - y_min) if y_max != y_min else 1.0
    ax.set_xlim(x_min - pad_x, x_max + pad_x)
    ax.set_ylim(y_min - pad_y, y_max + pad_y)

    plt.tight_layout()
    plt.show()



def plot_1d_loaded_structure(
    fem_nodes, fem_elements, u, stress_values, stress_type, ndof_per_node, scale=1.0, title=None, show_labels=False
):
    """
    Plots a 1D structure (beam or bar) with deformed shape and color-coded stress.
    """

    num_nodes = fem_nodes.shape[0]

    x = fem_nodes[:, 0]
    y_undeformed = np.zeros_like(x)

    # Extract displacements
    if ndof_per_node == 1:
        # Axial deformation
        x_deformed = x + scale * u
        y_deformed = y_undeformed
    elif ndof_per_node >= 2:
        # Vertical deformation (assume first DOF is vertical)
        u_y = np.array([u[i * ndof_per_node] for i in range(num_nodes)])
        x_deformed = x
        y_deformed = y_undeformed + scale * u_y
    else:
        raise ValueError("Unsupported number of DOFs per node")

    # Line segments and stress
    lines_undeformed = []
    lines_deformed = []
    stress_colors = []

    for e, (n1, n2) in enumerate(fem_elements):
        lines_undeformed.append([[x[n1], y_undeformed[n1]], [x[n2], y_undeformed[n2]]])
        lines_deformed.append([[x_deformed[n1], y_deformed[n1]], [x_deformed[n2], y_deformed[n2]]])
        stress_colors.append(stress_values[e])

    stress_colors = np.array(stress_colors)
    cmap = colormaps['viridis']
    color_norm = (stress_colors - stress_colors.min()) / (stress_colors.max() - stress_colors.min() + 1e-9)
    color_mapped = cmap(color_norm)

    # Plot
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.set_aspect('auto')

    ax.add_collection(LineCollection(lines_undeformed, linestyle='--',colors='gray', linewidths=2, label='Undeformed'))
    ax.add_collection(LineCollection(lines_deformed, colors=color_mapped, linewidths=3, label='Deformed'))

    # Optional stress labels
    if show_labels:
        for e, (n1, n2) in enumerate(fem_elements):
            xm = (x_deformed[n1] + x_deformed[n2]) / 2
            ym = (y_deformed[n1] + y_deformed[n2]) / 2
            ax.text(xm, ym, f"{stress_colors[e]:.1f}", fontsize=8, ha='center',
                    bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap)
    sm.set_array(stress_colors)
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label(f"{stress_type.capitalize()}[MPa]")

    # Set title and labels
    ax.set_title(title or f"1D Deformation and {stress_type.capitalize()}")
    ax.set_xlabel("X [mm]")
    ax.set_ylabel("Displacement [mm]")
    ax.grid(True)
    ax.legend()

    # Manual axis limits
    all_x = np.concatenate((x, x_deformed))
    all_y = np.concatenate((y_undeformed, y_deformed))
    x_pad = 0.05 * (all_x.max() - all_x.min()) or 1.0
    y_pad = 0.05 * (all_y.max() - all_y.min()) or 1.0

    ax.set_xlim(all_x.min() - x_pad, all_x.max() + x_pad)
    ax.set_ylim(all_y.min() - y_pad, all_y.max() + y_pad)

    plt.tight_layout()
    plt.show()






from matplotlib.animation import FuncAnimation
from IPython.display import HTML

def animate_mode_shape_2D(mode_index, eigenvecs, node_coords, free_dofs, K, elements, ndof_per_node, amplification=100, save_as=None):
    """
    Animate a mode shape for a 2D frame (with 3 DOFs per node).

    Parameters:
        mode_index: Index of the mode (0-based)
        eigenvecs: Array of eigenvectors from reduced system
        node_coords: (N, 2) array of node positions
        free_dofs: List of free DOF indices
        K: Global stiffness matrix (to get total DOFs)
        elements: (E, 2) array with node indices for each element
        amplification: Visual scaling factor for displacements
        save_as: Optional filename to save animation (e.g., "mode1.gif")
    
    Returns:
        HTML animation (for Jupyter notebooks)
    """
    num_nodes = node_coords.shape[0]
    total_dofs = K.shape[0]

    # Reconstruct full mode shape including constrained DOFs
    full_mode = np.zeros(total_dofs)
    full_mode[free_dofs] = eigenvecs[:, mode_index]

    # Split displacements into per-node ux and uy
    u_disp = full_mode[0::ndof_per_node]  # ux
    v_disp = full_mode[1::ndof_per_node]  # uy
    #theta  = full_mode[2::ndof_per_node]  # unused in plot, could rotate if desired

    # Normalize mode shape
    max_disp = np.max(np.sqrt(u_disp**2 + v_disp**2))
    u_disp /= max_disp
    v_disp /= max_disp

    u_disp *= amplification
    v_disp *= amplification

    # Undeformed shape
    x_static = node_coords[:, 0]
    y_static = node_coords[:, 1]

    # Setup plot
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect('equal')
    ax.set_title(f"Mode Shape {mode_index + 1}")
    ax.set_xlabel("X [mm]")
    ax.set_ylabel("Y [mm]")
    ax.grid(True)

    # Plot undeformed structure
    for el in elements:
        x = [node_coords[el[0], 0], node_coords[el[1], 0]]
        y = [node_coords[el[0], 1], node_coords[el[1], 1]]
        ax.plot(x, y, 'k--', linewidth=1, alpha=0.3)

    line_collection = []
    for _ in elements:
        line, = ax.plot([], [], 'b-', lw=2)
        line_collection.append(line)

    margin = amplification * 1.5
    ax.set_xlim(x_static.min() - margin, x_static.max() + margin)
    ax.set_ylim(y_static.min() - margin, y_static.max() + margin)

    # Animation update
    def update(frame):
        factor = np.sin(2 * np.pi * frame / 60)
        displaced_coords = node_coords + factor * np.column_stack((u_disp, v_disp))

        for i, el in enumerate(elements):
            x = [displaced_coords[el[0], 0], displaced_coords[el[1], 0]]
            y = [displaced_coords[el[0], 1], displaced_coords[el[1], 1]]
            line_collection[i].set_data(x, y)

        return line_collection

    ani = FuncAnimation(fig, update, frames=60, interval=50, blit=True)

    if save_as:
        ani.save(save_as, writer="pillow")
        print(f"Animation saved as {save_as}")
        plt.close(fig)
    else:
        plt.close(fig)
        return HTML(ani.to_jshtml())


    
def animate_mode_shape_1D(mode_index, eigenvecs, node_coords, free_dofs, K, elements, ndof_per_node, amplification=100, save_as=None):
    """
    Animate a mode shape for 1D bar or beam elements with moving node labels.

    Parameters:
        mode_index: Index of the mode (0-based)
        eigenvecs: Array of eigenvectors from reduced system
        node_coords: (N, 2) array of node positions (x, y)
        free_dofs: List of free DOF indices
        K: Global stiffness matrix (used to get total DOFs)
        elements: (E, 2) array with node indices for each element
        ndof_per_node: Number of DOFs per node (1 for bar, 2 for beam)
        amplification: Scaling factor for visual displacement
        save_as: Optional filename to save animation (e.g., "mode1.gif")

    Returns:
        HTML animation (for Jupyter notebooks) if not saving
    """
    total_dofs = K.shape[0]
    num_nodes = node_coords.shape[0]

    full_mode = np.zeros(total_dofs)
    full_mode[free_dofs] = eigenvecs[:, mode_index]

    if ndof_per_node == 1:
        u_disp = full_mode[0::ndof_per_node]
        v_disp = np.zeros_like(u_disp)
    elif ndof_per_node == 2:
        u_disp = np.zeros(num_nodes)
        v_disp = full_mode[0::ndof_per_node]
    else:
        raise ValueError("Only 1 or 2 DOFs per node supported for 1D elements.")

    max_disp = np.max(np.sqrt(u_disp**2 + v_disp**2))
    if max_disp == 0:
        max_disp = 1.0
    u_disp /= max_disp
    v_disp /= max_disp

    u_disp *= amplification
    v_disp *= amplification

    x_static = node_coords[:, 0]
    y_static = node_coords[:, 1]

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.set_aspect('equal')
    ax.set_title(f"Mode Shape {mode_index + 1}")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.grid(True)

    # Plot undeformed structure (gray dashed lines)
    for el in elements:
        x = [node_coords[el[0], 0], node_coords[el[1], 0]]
        y = [node_coords[el[0], 1], node_coords[el[1], 1]]
        ax.plot(x, y, 'k--', linewidth=1, alpha=0.3)

    # Initialize animated lines
    line_collection = []
    for _ in elements:
        line, = ax.plot([], [], 'b-', lw=2)
        line_collection.append(line)

    # Initialize moving nodes and text labels
    node_dots = ax.plot([], [], 'ro', markersize=4)[0]
    node_texts = [ax.text(0, 0, str(i), color='red', fontsize=8, ha='left', va='bottom') for i in range(num_nodes)]

    # Set limits
    margin_x = amplification * 1.5
    margin_y = amplification * 1.5
    ax.set_xlim(x_static.min() - margin_x, x_static.max() + margin_x)
    ax.set_ylim(y_static.min() - margin_y, y_static.max() + margin_y)

    def update(frame):
        factor = np.sin(2 * np.pi * frame / 60)
        displaced_coords = node_coords + factor * np.column_stack((u_disp, v_disp))

        # Update lines
        for i, el in enumerate(elements):
            x = [displaced_coords[el[0], 0], displaced_coords[el[1], 0]]
            y = [displaced_coords[el[0], 1], displaced_coords[el[1], 1]]
            line_collection[i].set_data(x, y)

        # Update node dots
        node_dots.set_data(displaced_coords[:, 0], displaced_coords[:, 1])

        # Update node text labels
        for i, txt in enumerate(node_texts):
            txt.set_position((displaced_coords[i, 0] + 0.02 * amplification, 
                              displaced_coords[i, 1] + 0.02 * amplification))

        return line_collection + [node_dots] + node_texts

    ani = FuncAnimation(fig, update, frames=60, interval=50, blit=True)

    if save_as:
        ani.save(save_as, writer="pillow")
        print(f"Animation saved as {save_as}")
        plt.close(fig)
    else:
        plt.close(fig)
        return HTML(ani.to_jshtml())
