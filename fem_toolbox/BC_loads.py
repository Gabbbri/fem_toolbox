
import re   # for parsing, regular expressions
import numpy as np
import matplotlib.pyplot as plt

def read_bc_and_forces(file_path):
    """
    Reads boundary conditions and forces from a structured file with compact BC syntax.

    Returns:
        bc_nodes: list of nodes with BCs
        bc_dofs: corresponding DOFs
        bc_values: corresponding values
        f_nodes: list of nodes with forces
        f_dofs: corresponding DOFs
        f_values: corresponding values
    """
    bc_nodes = []
    bc_dofs = []
    bc_values = []

    f_nodes = []
    f_dofs = []
    f_values = []

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            if line.startswith("BC"):

                # extract node number and (dof value) pairs
                node_match = re.match(r"BC\s+(\d+)", line)
                if not node_match:
                    raise ValueError(f"Invalid BC line format: {line}")
                node = int(node_match.group(1))

                # find all (dof value) pairs
                pairs = re.findall(r"\((\d+)\s+([-+]?\d*\.?\d+)\)", line)
                for dof_str, val_str in pairs:
                    dof = int(dof_str)
                    value = float(val_str)

                    bc_nodes.append(node)
                    bc_dofs.append(dof)
                    bc_values.append(value)

            elif line.startswith("FORCE"):
                parts = line.split()
                if len(parts) != 4:
                    raise ValueError(f"Invalid FORCE line: {line}")
                node = int(parts[1])
                dof = int(parts[2])
                value = float(parts[3])

                f_nodes.append(node)
                f_dofs.append(dof)
                f_values.append(value)

            else:
                raise ValueError(f"Unrecognized line: {line}")
            
    # get unique constrained nodes
    constrained_nodes = sorted(set(bc_nodes))

    return (
        constrained_nodes,
        bc_nodes,
        bc_dofs,
        bc_values,
        f_nodes,
        f_dofs,
        f_values
    )



def validate_constraints(num_nodes, num_dofs, bc_nodes, bc_dofs):
    """
    Validates that the structure is not under-constrained.
    Validates that the nodes for external loads and BC exist

    Args:
        num_nodes (int): total number of nodes
        num_dofs (int) : num of dofs per node (bar -> 1dof, 2d beam -> 3 dofs)
        bc_nodes (list of int): list of nodes with BCs (can be repeated)
        bc_dofs (list of int): corresponding DOFs (0=u, 1=v, 2=theta)

    Raises:
        ValueError: if the structure is not adequately constrained
    """
    total_dofs = num_nodes * num_dofs
    constrained_dof_ids = [node * num_dofs + dof for node, dof in zip(bc_nodes, bc_dofs)]
    constrained_dof_ids = sorted(set(constrained_dof_ids))  # unique constrained dof ids

    for dof in constrained_dof_ids:
        if dof > total_dofs:
            raise ValueError("Trying to constrain dofs that don't exist! Check geometry and boundary conditions file")

    num_constrained = len(constrained_dof_ids)

    automatically_constrained = 0

    if num_dofs == 3 and  num_constrained < 3:
        raise ValueError("Structure is under-constrained: fewer than 3 total constrained DOFs.")

    if num_dofs == 2:
        #automatically_constrained = 1
        if num_constrained < 2:
            raise ValueError("Structure is under-constrained: fewer than 3 total constrained DOFs.")

    if num_dofs == 1 and any(i >= 1 for i in bc_dofs):  # 2D structure with bar elements
        #automatically_constrained = 2
        if num_constrained < 2:
            raise ValueError("Structure is under-constrained: fewer than 2 total constrained DOFs.")
        
    elif num_dofs == 1:  # 1D structure with bar elements
        #automatically_constrained = 2
        if num_constrained < 1:
            raise ValueError("Structure is under-constrained: fewer than 1 total constrained DOFs.")

    if num_constrained >= total_dofs:
        raise ValueError("Structure is over-constrained: all DOFs are fixed (structure cannot deform).")

    print(f"Constraint check passed: {num_constrained} DOFs constrained out of {total_dofs}.")




def plot_fem_model(node_coords, fem_elements,bc_nodes, bc_dofs, bc_vals, force_nodes, force_dofs, force_vals, scale_force=0.1, scale_moment=0.1):
    """
    Visualizes the FEM structure, boundary conditions, and forces.
    - node_coords: (N, 2) array of node positions
    - fem_elements: (E, 2) array of element node indices
    - *_nodes: lists of node indices
    - *_dofs: lists of DOFs (0=u, 1=v, 2=theta)
    - *_vals: values of BCs or forces
    """
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Plot elements
    for n1, n2 in fem_elements:
        x = [node_coords[n1, 0], node_coords[n2, 0]]
        y = [node_coords[n1, 1], node_coords[n2, 1]]
        ax.plot(x, y, 'k-', lw=2)

    # Plot nodes with labels
    for idx, (x, y) in enumerate(node_coords):
        ax.plot(x, y, 'ro', markersize=4)
        ax.text(x + 0.02, y + 0.02, str(idx), fontsize=9, color='blue')

    # Plot boundary conditions
    for node, dof, val in zip(bc_nodes, bc_dofs, bc_vals):
        x, y = node_coords[node]
        if dof == 0:  # u
            ax.plot(x - 50*scale_force, y, 'bs', markersize=10, label='u=0' if val == 0 else None)
        elif dof == 1:  # v
            ax.plot(x, y - 50*scale_force, 'gs', markersize=10, label='v=0' if val == 0 else None)
        elif dof == 2:  # theta
            ax.plot(x, y, 'ms', markersize=10, label='θ=0' if val == 0 else None)

    # Plot forces
    for node, dof, val in zip(force_nodes, force_dofs, force_vals):
        x, y = node_coords[node]
        if dof == 0:  # u
            dx = val*scale_force
            ax.arrow(x, y, dx, 0, head_width=0.2*np.abs(val)*scale_force, color='blue', length_includes_head=True)
        elif dof == 1:  # v
            dy = val*scale_force
            ax.arrow(x, y, 0, dy, head_width=0.2*np.abs(val)*scale_force, color='green', length_includes_head=True)
        elif dof == 2:  # theta (moment)
            radius = scale_moment
            theta = np.linspace(0, 2 * np.pi, 100)
            ax.plot(x + radius * np.cos(theta), y + radius * np.sin(theta), 'r--')
            ax.text(x + radius, y + radius, 'M', fontsize=9, color='red')

    ax.set_aspect('equal')
    ax.set_title("FEM Geometry with BCs and Forces")
    ax.grid(True)
    ax.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.tight_layout()
    plt.show()
