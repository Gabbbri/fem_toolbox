
import numpy as np
import matplotlib.pyplot as plt

def discretize(node_coords_beam, beam_connectivity, beam_crossSections, elements_per_beam=5, strategy="uniform"):
    """
    Discretizes beams into FEM elements, keeping A and I per element.
    
    Returns:
        fem_nodes : (N, 2) ndarray
        fem_elements : list of (node_i, node_j)
        element_crossSection : list of dicts with A and I for each element
    """
    fem_nodes = []
    node_index_map = {}  # maps (x,y) to node index ((x, y) = key, index = value)
    fem_elements = []
    element_crossSection = []

    # avoid adding duplicates to the node list, and ensures to get consistent indexing
    current_index = 0
    def get_or_add_node(x, y):
        nonlocal current_index
        key = (round(x, 5), round(y, 5))
        if key not in node_index_map:
            node_index_map[key] = current_index
            fem_nodes.append([x, y])
            current_index += 1
        return node_index_map[key]

    
    for beam_idx, (n1, n2) in enumerate(beam_connectivity):
        p1 = node_coords_beam[n1]
        p2 = node_coords_beam[n2]

        A = beam_crossSections[beam_idx]["A"]
        I = beam_crossSections[beam_idx]["I"]


        # uniform subdivision of the beam into elements
        for i in range(elements_per_beam):
            # t1,2: parametric values, "how far we are along the beam"
            t1 = i / elements_per_beam
            t2 = (i + 1) / elements_per_beam

            # linear interpolation along the element
            x1, y1 = (1 - t1) * p1 + t1 * p2
            x2, y2 = (1 - t2) * p1 + t2 * p2

            idx1 = get_or_add_node(x1, y1)
            idx2 = get_or_add_node(x2, y2)

            fem_elements.append([idx1, idx2])
            element_crossSection.append({"A": A, "I": I})

    return np.array(fem_nodes), np.array(fem_elements), element_crossSection 



def plot_discretized_geometry_2D(original_nodes, original_beams, fem_nodes, fem_elements, show_elements=True, show_nodes=True, node_label_offset=0.02):
    """
    Visualizes the original beam geometry and the discretized FEM mesh.

    Parameters:
        original_nodes : (N, 2) array of unique input nodes
        original_beams : (B, 2) array of beam connectivity (node indices)
        fem_nodes : (M, 2) array of discretized FEM nodes
        fem_elements : (E, 2) array of FEM element connectivity
        show_elements : bool - whether to show element numbers
        show_nodes : bool - whether to show node numbers
        node_label_offset : float - visual offset for node label
    """
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Plot original geometry (beam lines)
    for start, end in original_beams:
        x_vals = [original_nodes[start][0], original_nodes[end][0]]
        y_vals = [original_nodes[start][1], original_nodes[end][1]]
        ax.plot(x_vals, y_vals, 'k--', linewidth=1, alpha=0.3)

    # Plot discretized FEM elements
    for i, (n1, n2) in enumerate(fem_elements):
        x_vals = [fem_nodes[n1][0], fem_nodes[n2][0]]
        y_vals = [fem_nodes[n1][1], fem_nodes[n2][1]]
        ax.plot(x_vals, y_vals, 'b-', linewidth=2)

        # Optional element ID at midpoint
        if show_elements:
            xm, ym = (np.array(x_vals).mean(), np.array(y_vals).mean())
            ax.text(xm, ym, str(i), color='blue', fontsize=8, ha='center', va='center')

    # Plot and label nodes
    for i, (x, y) in enumerate(fem_nodes):
        ax.plot(x, y, 'ro', markersize=4)
        if show_nodes:
            ax.text(x + node_label_offset, y + node_label_offset, str(i), color='red', fontsize=9)

    ax.set_aspect('equal')
    ax.set_title("Discretized Geometry")
    ax.grid(True)
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.show()



def plot_discretized_geometry_1D(original_nodes, original_beams, fem_nodes, fem_elements, show_elements=True, show_nodes=True, node_label_offset=0.02):
    """
    Visualizes the original 1D beam geometry and the discretized FEM mesh for 1D structures.
    
    Parameters:
        original_nodes : (N, 2) array of unique input nodes
        original_beams : (B, 2) array of beam connectivity (node indices)
        fem_nodes : (M, 2) array of discretized FEM nodes
        fem_elements : (E, 2) array of FEM element connectivity
        show_elements : bool - whether to show element numbers
        show_nodes : bool - whether to show node numbers
        node_label_offset : float - visual offset for node label
    """
    fig, ax = plt.subplots(figsize=(8, 2))  # We keep the height smaller for 1D geometry
    
    # Plot original geometry (beam lines), only in the horizontal direction
    for start, end in original_beams:
        x_vals = [original_nodes[start][0], original_nodes[end][0]]
        y_vals = [original_nodes[start][1], original_nodes[end][1]]
        ax.plot(x_vals, y_vals, 'k--', linewidth=1, alpha=0.3)
    
    # Plot discretized FEM elements, only in the horizontal direction
    for i, (n1, n2) in enumerate(fem_elements):
        x_vals = [fem_nodes[n1][0], fem_nodes[n2][0]]  # Horizontal displacement (x-axis)
        y_vals = [fem_nodes[n1][1], fem_nodes[n2][1]]  # Vertical displacement (should be 0 for 1D)
        
        # Only plot the horizontal components
        ax.plot(x_vals, y_vals, 'b-', linewidth=2)

        # Optional element ID at midpoint
        if show_elements:
            xm, ym = (np.array(x_vals).mean(), np.array(y_vals).mean())
            ax.text(xm, ym, str(i), color='blue', fontsize=8, ha='center', va='center')

    # Plot and label nodes
    for i, (x, y) in enumerate(fem_nodes):
        ax.plot(x, y, 'ro', markersize=4)
        if show_nodes:
            ax.text(x + node_label_offset, y + node_label_offset, str(i), color='red', fontsize=9)

    # Set aspect for 1D structure
    ax.set_aspect('auto')
    ax.set_title("Discretized Geometry (1D)")
    ax.grid(True)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    plt.show()

