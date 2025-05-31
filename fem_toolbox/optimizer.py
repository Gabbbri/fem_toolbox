from scipy.optimize import minimize
import numpy as np
import matplotlib.pyplot as plt
import fem_toolbox as ft

def optimize_crossSections4stress_trusses(
    element_crossSections,
    sigma_max,
    fem_nodes,
    fem_elements,
    mat_properties,
    ndof_per_node,
    bc_nodes,
    bc_dofs,
    bc_values,
    f_nodes,
    f_dofs,
    f_values,
    assembleK,
    assembleF,
    k_local_func,
    rotation_func,
    A_min=10,   # units are mm
    A_max=1e5,
    section_shape="rectangle"
):
    num_elements = len(element_crossSections)
    original_inertias = [e['I'] for e in element_crossSections]

    def objective(A):   # wight minimization (= volume minimization for constant density)
        lengths = []
        for i, j in fem_elements:
            xi, yi = fem_nodes[i]
            xj, yj = fem_nodes[j]
            L = np.sqrt((xj - xi)**2 + (yj - yi)**2)
            lengths.append(L)
        lengths = np.array(lengths)

        element_volume = A * lengths  # Element-wise multiplication
        return np.sum(element_volume)


    #def objective(A):
    #    return np.sum(A)

    def constraint_fun(A):
        # Convert back A vector to list of dictionaries (all fem_toolbox functions work with that)
        element_areas = [{'A': A[i], 'I': original_inertias[i]} for i in range(num_elements)]

        # Assemble stiffness matrix and force vector
        K = assembleK(k_local_func, rotation_func, fem_nodes, fem_elements,
                      element_areas, mat_properties, ndof_per_node)
        total_dofs = K.shape[0]

        struct2D = False
        if rotation_func is not None:
            struct2D = True

        f = assembleF(f_nodes, f_dofs, f_values, total_dofs, ndof_per_node, struct2D_trussElements=struct2D)

        # Solve static system
        u = ft.femsolver.static_analysis(K, f, bc_nodes, bc_dofs, bc_values, ndof_per_node)

        # Evaluate stresses
        _, axial_stress, _, _, _, _ = ft.postprocessing.eval_stress(
            k_local_func, rotation_func, u, fem_elements, fem_nodes,
            element_areas, mat_properties, section_shape, ndof_per_node
        )

        # Constraint: sigma_i <= sigma_max ⇒ sigma_max - sigma_i >= 0
        return sigma_max - np.abs(axial_stress)

    # Define bounds and constraints for SLSQP
    bounds = [(A_min, A_max) for _ in range(num_elements)]
    constraints = {
        'type': 'ineq',
        'fun': constraint_fun  # if I don't give jacobian, SciPy estimates it via finite differences
    }

    x0 = np.array([e['A'] for e in element_crossSections])

    result = minimize(      # using built-in finite differences
        fun=objective,
        x0=x0,
        method='SLSQP',
        bounds=bounds,
        constraints=[constraints],
        options={'disp': True, 'maxiter': 100}
    )

    return result.x


def checkOptimization_stresses_trusses(
    element_crossSections,
    optimized_areas,
    max_stress,
    fem_nodes,
    fem_elements,
    mat_properties,
    ndof_per_node,
    bc_nodes,
    bc_dofs,
    bc_values,
    f_ext,
    f_nodes,
    f_dofs,
    assembleK,
    assembleF,
    k_local_func,
    rotation_func,
    section_shape="rectangle",
    verbose=True,
):
    """
    Compute and compare the axial stresses in the structure before and after optimization.

    Parameters
    ----------
    original_areas : list of dictionaries, each with keys 'A' and 'I'.
    optimized_areas : list of dictionaries, same structure.
    max_stress : float, maximum allowable axial stress.
    ... : FEM assembly parameters.
    verbose : bool, whether to print formatted output.
    """
    num_elements = len(fem_elements)
    
    # building the usual list of dictionaries for the area properties
    original_inertias = [e['I'] for e in element_crossSections]
    optimized_areas = [{'A': A, 'I': I} for A, I in zip(optimized_areas, original_inertias)]

    def compute_stresses(areas):
        # Assemble global stiffness matrix
        K = assembleK(k_local_func, rotation_func, fem_nodes, fem_elements,
                      areas, mat_properties, ndof_per_node)
        total_dofs = K.shape[0]

        struct2D = False
        if rotation_func is not None:
            struct2D = True
        # Assemble load vector
        F = assembleF(f_nodes, f_dofs, f_ext, total_dofs, ndof_per_node, struct2D_trussElements=False)

        # Solve static system
        u = ft.femsolver.static_analysis(K, f_ext, bc_nodes, bc_dofs, bc_values, ndof_per_node)

        # Evaluate stresses
        _, axial_stress, _, _, _, _ = ft.postprocessing.eval_stress(
            k_local_func, rotation_func, u, fem_elements, fem_nodes,
            areas, mat_properties, section_shape, ndof_per_node
        )

        return axial_stress

    # Compute
    orig_stresses = compute_stresses(element_crossSections)
    opt_stresses = compute_stresses(optimized_areas)

    if verbose:
        print("Element  | Area (orig) | Area (opt) | Stress (orig) | Stress (opt) | max_stress")
        print("---------|-------------|------------|----------------|----------------|----------")
        for i in range(num_elements):
            A_orig = element_crossSections[i]['A']
            A_opt = optimized_areas[i]['A']
            print(f"{i:^8} | {A_orig:^11.2f} | {A_opt:^10.2f} |"
                  f" {orig_stresses[i]:^14.2f} | {opt_stresses[i]:^14.2f} | {max_stress:>9.2f}")

    #return orig_stresses, opt_stresses



def optimize_crossSections4frequency_trusses(
    element_crossSections,
    forbidden_range,  # tuple (freq_min, freq_max)
    fem_nodes,
    fem_elements,
    mat_properties,
    ndof_per_node,
    bc_nodes,
    bc_dofs,
    bc_values,
    assembleK,
    assembleM,
    num_modes,
    k_local_func,
    m_local_func,
    rotation_func=None,
    A_min=50,
    A_max=3000,
    margin=0.5,
    uniformity_penalty_weight=1e-1,  # Set > 0 to enable uniformity penalty
    smoothness_penalty_weight=0.3e-1,
    boundary_penalty_weight=100,
    initial_perturbation=0.2        # Set > 0 to add initial noise (e.g., 0.05 for ±5%)
):


    num_elements = len(element_crossSections)
    original_inertias = [e['I'] for e in element_crossSections]
    freq_min, freq_max = forbidden_range
    margin = 1e-2  # safety margin



    # --- Initial area vector, PERTURBED ---
    x0 = np.array([e['A'] for e in element_crossSections])

    if np.average(x0) < A_min:
        x0 = np.full_like(x0, A_min + np.average(x0))
    if np.average(x0) >= A_max:
        x0 = np.full_like(x0, A_max - np.average(x0) )
    x0 = np.clip(x0, A_min, A_max)
    if initial_perturbation > 0.0:
        perturbation = np.random.uniform(-initial_perturbation, initial_perturbation, size=len(x0))
        x0 *= (1.0 + perturbation)
            
        x0 = np.clip(x0, A_min, A_max)
        print(x0)

    # --- Bounds ---
    bounds = [(A_min, A_max) for _ in x0]


    def objective(A):   # wight minimization (= volume minimization for constant density)
       # lengths = []
       # for i, j in fem_elements:
       #     xi, yi = fem_nodes[i]
       #     xj, yj = fem_nodes[j]
       #     L = np.sqrt((xj - xi)**2 + (yj - yi)**2)
       #     lengths.append(L)
#
       # lengths = np.array(lengths)
       # element_volume = A * lengths  # Element-wise multiplication
       # structure_volume = np.sum(element_volume)
        structure_volume = np.sum(A)
#
        if uniformity_penalty_weight > 0.0:
            global_variance = np.var(A)
            local_jump_penalty = np.sum(np.diff(A)**2)

            bound_penalty = np.sum((A - A_min) < 1) * boundary_penalty_weight  # Penalize sticking to min bound

            penalty = (-uniformity_penalty_weight * global_variance +
                   smoothness_penalty_weight * local_jump_penalty +
                   bound_penalty)

            return structure_volume + penalty

        return structure_volume

    def constraint_fun(A):
        # Convert A vector to list of dictionaries
        element_areas = [{'A': A[i], 'I': original_inertias[i]} for i in range(num_elements)]

        # Assemble K and M
        K = assembleK(k_local_func, rotation_func, fem_nodes, fem_elements,
                      element_areas, mat_properties, ndof_per_node)
        M = assembleM(m_local_func, rotation_func, fem_nodes, fem_elements,
                      element_areas, mat_properties, ndof_per_node)

        # Solve the frequencies
        frequencies, _ , _ = ft.femsolver.modal_analysis(K, M, bc_nodes, bc_dofs, ndof_per_node, num_modes = num_modes, verbose=False)

        # Impose constraints
        constraints = []

        for frequency in frequencies:
            if freq_min + margin < frequency < freq_max - margin:
                # Inside the forbidden band -> violation
                # calculate distances from frequencies to domain boundaries. You impose them < 0. Since scipy inequalities want > 0, I set  (- distances) > 0 
                constraints.append(-(frequency - freq_min - margin))  # want <= 0
                constraints.append(-(freq_max - frequency - margin))  # want <= 0
            else:
                # Outside forbidden band → constraint already satisfied
                constraints.append(1.0)  # dummy positive value
                constraints.append(1.0)

        return np.array(constraints)
    

    constraints = {
        'type': 'ineq',
        'fun': constraint_fun  # Finite-diff sensitivities used automatically
    }

    # --- Optimization ---
    result = minimize(
        fun=objective,
        x0=x0,
        method='SLSQP',
        #method='trust-constr',     # definitely to not use
        bounds=bounds,
        constraints=[constraints],
        options={'disp': True, 'maxiter': 1000}
    )

    if not result.success:
        print("\n\t\
              Optimization was not successfull. Please try again\n")

    return result, result.x


import random

def optimizeFrequency_random_hyperparameter_search(
    max_trials,
    element_crossSections,
    forbidden_range,
    fem_nodes,
    fem_elements,
    mat_properties,
    ndof_per_node,
    bc_nodes,
    bc_dofs,
    bc_values,
    assembleK,
    assembleM,
    num_modes,
    k_local_func,
    m_local_func,
    A_min,
    A_max,
    rotation_func
):
    for trial in range(max_trials):
        # Random hyperparameters in sensible ranges
        uniformity_penalty = random.uniform(0.2, 0.4)
        smoothness_penalty = random.uniform(0.1, 0.2)
        boundary_penalty = random.uniform(100, 400.0)

        print(f"\n[Trial {trial+1}] Trying with:")
        print(f"  Uniformity penalty: {uniformity_penalty:.4f}")
        print(f"  Smoothness penalty: {smoothness_penalty:.4f}")
        print(f"  Boundary penalty:   {boundary_penalty:.2f}")

        result, _ = ft.optimizer.optimize_crossSections4frequency_trusses(
            element_crossSections=element_crossSections,
            forbidden_range=forbidden_range,
            fem_nodes=fem_nodes,
            fem_elements=fem_elements,
            mat_properties=mat_properties,
            ndof_per_node=ndof_per_node,
            bc_nodes=bc_nodes,
            bc_dofs=bc_dofs,
            bc_values=bc_values,
            assembleK=assembleK,
            assembleM=assembleM,
            num_modes=num_modes,
            k_local_func=k_local_func,
            m_local_func=m_local_func,
            rotation_func=rotation_func,
            A_min=A_min,
            A_max=A_max,
            uniformity_penalty_weight=uniformity_penalty,
            smoothness_penalty_weight=smoothness_penalty,
            boundary_penalty_weight=boundary_penalty,
        )

        if result.success:
            print("\n✅ Success! Exiting hyperparameter tuning.")
            return result.x, {
                'uniformity_penalty_weight': uniformity_penalty,
                'smoothness_penalty_weight': smoothness_penalty,
                'boundary_penalty_weight': boundary_penalty
            }

        print(f"❌ Trial {trial+1} failed: {result.message}")

    print("\n⚠️ No successful optimization found in max_trials.")
    
    return None, None


def checkOptimization_frequencies_trusses(
    element_crossSections,
    optimized_area_values,
    forbidden_frequency_band,
    fem_nodes,
    fem_elements,
    mat_properties,
    ndof_per_node,
    bc_nodes,
    bc_dofs,
    assembleK,
    assembleM,
    k_local_func,
    m_local_func,
    rotation_func=None,
    num_modes=5,
    verbose=True,
):
    """
    Compare modal frequencies before and after optimization.

    Parameters
    ----------
    original_areas : list of dicts, each with 'A' and 'I'
    optimized_areas : list of dicts, same structure
    num_modes : int, number of modes to compute
    verbose : bool, whether to print formatted output
    """

    # building the usual list of dictionaries for the area properties
    original_inertias = [e['I'] for e in element_crossSections]
    optimized_area_values = [{'A': A_val, 'I': I} for A_val, I in zip(optimized_area_values, original_inertias)]

    def compute_frequencies(area_data):
        K = assembleK(k_local_func, rotation_func, fem_nodes, fem_elements,
                      area_data, mat_properties, ndof_per_node)
        M = assembleM(m_local_func, rotation_func, fem_nodes, fem_elements,
                      area_data, mat_properties, ndof_per_node)
        freqs, _, _ = ft.femsolver.modal_analysis(K, M, bc_nodes, bc_dofs, ndof_per_node, num_modes=num_modes, verbose=False)
        return freqs

    freqs_orig = compute_frequencies(element_crossSections)
    freqs_opt = compute_frequencies(optimized_area_values)

    if verbose:
        print("\nElement-wise Area Comparison:")
        print("Element | Area (original) | Area (optimized)")
        print("--------|------------------|------------------")
        for i, (a1, a2) in enumerate(zip(element_crossSections, optimized_area_values)):
            print(f"{i:^8} | {a1['A']:^16.2f} | {a2['A']:^16.2f}")

        print(f"\nForbidden frequency band\n    f_min = {forbidden_frequency_band[0]}\n    f_max = {forbidden_frequency_band[1]}\n")

        print("\nModal Frequencies Comparison:")
        print("Mode    | Frequency (orig) | Frequency (opt)")
        print("--------|------------------|-----------------")
        for i, (f1, f2) in enumerate(zip(freqs_orig, freqs_opt), 1):
            print(f"{i:^8} | {f1:^16.4f} | {f2:^15.4f}")

    #return freqs_orig, freqs_opt

def optim_plot_1D_beam_geometry(
    fem_nodes, fem_elements, areas, section_shape="rectangle", color='skyblue',
    padding_factor=1.2, min_thickness=1e-2, figsize=(10, 6), show_labels=True
):
    """
    Plot a visual representation of 1D beam geometry based on cross-sectional areas.

    Arguments:
        fem_nodes : (N, 2) array of node coordinates
        fem_elements : list of (node_i, node_j) tuples
        areas : list or array of cross-sectional areas
        section_shape : "rectangle" or "circular"
        color : color of the filled cross-section bars
        padding_factor : extra space added around the plot
        min_thickness : minimum visual thickness to avoid collapsing elements
        figsize : figure size in inches (width, height)
        show_labels : whether to annotate the beam with thickness/diameter labels
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    all_x = []
    all_y = []

    for idx, (i, j) in enumerate(fem_elements):
        x0, y0 = fem_nodes[i]
        x1, y1 = fem_nodes[j]
        A = areas[idx]
        
        # Convert area to visual width
        if section_shape == "rectangle":
            width = max(np.sqrt(A), min_thickness)
            label_text = f"{width:.1f} mm"
        elif section_shape == "circular":
            diameter = max(2 * np.sqrt(A / np.pi), min_thickness)
            width = diameter
            label_text = f"Ø {diameter:.1f} mm"
        else:
            raise ValueError("Unsupported section shape")

        # Direction and normal vector
        dx, dy = x1 - x0, y1 - y0
        length = np.hypot(dx, dy)
        if length == 0:
            continue
        nx, ny = -dy / length, dx / length  # normal vector

        # Rectangle polygon
        offset_x = 7 * width * nx
        offset_y = 7 * width * ny
        polygon = np.array([
            [x0 - offset_x, y0 - offset_y],
            [x0 + offset_x, y0 + offset_y],
            [x1 + offset_x, y1 + offset_y],
            [x1 - offset_x, y1 - offset_y]
        ])

        ax.fill(polygon[:, 0], polygon[:, 1], color=color, edgecolor='k', linewidth=0.8)
        all_x.extend(polygon[:, 0])
        all_y.extend(polygon[:, 1])

        # Midpoint for label
        if show_labels:
            xm = 0.5 * (x0 + x1)
            ym = 0.5 * (y0 + y1)
            ax.text(xm, ym, label_text, fontsize=8, ha='center', va='center', color='black', bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, boxstyle='round,pad=0.2'))

    if not all_x:
        ax.text(0.5, 0.5, "No valid elements to display", ha='center', va='center')
        return

    # Axis settings
    x_min, x_max = min(all_x), max(all_x)
    y_min, y_max = min(all_y), max(all_y)
    x_pad = (x_max - x_min) * (padding_factor - 1)
    y_pad = (y_max - y_min) * (padding_factor - 1)

    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel("X", fontsize=12)
    ax.set_ylabel("Y", fontsize=12)
    ax.set_title(f"Optimized Beam Geometry (Section: {section_shape})", fontsize=14)
    ax.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    ax.get_yaxis().set_visible(False)

    plt.show()


def optim_plot_truss_geometry(
    fem_nodes,
    fem_elements,
    areas,
    optim="stress",
    section_shape="rectangle",   # or "circular"
    color='steelblue',
    min_thickness=0.5,
    scale_factor=1.0,
    show_labels=True,
    figsize=(10, 8),
    padding_factor=1.1
):
    """
    Plot truss geometry with element thickness scaled by cross-sectional area (lines only).

    Arguments:
        fem_nodes : array of (x, y) coordinates
        fem_elements : list of (start_node_idx, end_node_idx) tuples
        areas : list of cross-sectional areas
        section_shape : "rectangle" or "circular"
        color : line color
        min_thickness : minimum visible linewidth
        max_linewidth : upper cap on linewidth (for visual consistency)
        scale_factor : visual scaling for line width
        show_labels : annotate elements with dimension labels
        figsize : size of the figure
        padding_factor : zoom padding multiplier
    """
    fig, ax = plt.subplots(figsize=figsize)
    fem_nodes = np.array(fem_nodes)

    all_x, all_y = [], []

    for idx, (i, j) in enumerate(fem_elements):
        x0, y0 = fem_nodes[i]
        x1, y1 = fem_nodes[j]
        A = areas[idx]

        # Compute physical dimension
        if section_shape == "rectangle":
            width = max(np.sqrt(A), min_thickness)
            label_text = f"{width:.1f} mm"
        elif section_shape == "circular":
            diameter = max(2 * np.sqrt(A / np.pi), min_thickness)
            width = diameter
            label_text = f"Ø {diameter:.1f} mm"
        else:
            raise ValueError("Unsupported section shape")

        # Ssale line thickness visually
        linewidth = width * scale_factor

        ax.plot([x0, x1], [y0, y1], color=color, linewidth=linewidth, solid_capstyle='round')
        all_x.extend([x0, x1])
        all_y.extend([y0, y1])

        # label 
        if show_labels:
            xm, ym = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
            ax.text(xm, ym, label_text, fontsize=8, ha='center', va='center',
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, boxstyle='round,pad=0.2'))

    if not all_x:
        ax.text(0.5, 0.5, "No valid elements to display", ha='center', va='center')
        return

    # put some axis limits (padding)
    x_min, x_max = min(all_x), max(all_x)
    y_min, y_max = min(all_y), max(all_y)
    x_pad = (x_max - x_min) * (padding_factor - 1)
    y_pad = (y_max - y_min) * (padding_factor - 1)

    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)
    ax.set_aspect('equal', adjustable='box')
    ax.set_title(f"Truss Geometry, optimized for {optim} constraint - Section: {section_shape}", fontsize=14)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()

