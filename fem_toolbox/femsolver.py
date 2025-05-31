import numpy as np
import fem_toolbox.elements
import inspect


def dof_map(ndofs_per_node, n1, n2):
    if ndofs_per_node == 1:     # bar element
        dof_map = [n1 * ndofs_per_node, n2 * ndofs_per_node]

    elif ndofs_per_node == 2:   # beam element
        dof_map = [
            n1 * ndofs_per_node,
            n1 * ndofs_per_node + 1,
            n2 * ndofs_per_node,
            n2 * ndofs_per_node + 1
        ]
        
    elif ndofs_per_node == 3:   # 2D beam/frame element
        dof_map = [
            ndofs_per_node * n1 + 0,
            ndofs_per_node * n1 + 1,
            ndofs_per_node * n1 + 2,
            ndofs_per_node * n2 + 0,
            ndofs_per_node * n2 + 1,
            ndofs_per_node * n2 + 2
        ]

    else:
        raise ValueError("Unsupported number of DOFs per node")

    return dof_map


def assembleK(k_local_func, rotation_func, fem_nodes, fem_elements, element_crossSections, mat_properties, ndof_per_node):
    """
    Assemble the global stiffness matrix.
    
    Arguments:
        k_local_func : pointer to the function for local stiffness matrix evaluation
        rotation_func : pointer to the function for local matrix rotation
        fem_nodes : (N, 2) array of node coordinates
        fem_elements : list of (node_i, node_j)
        element_crossSections : list of dicts with 'A' and 'I'
        element_materials : list of dicts with 'E' (and possibly 'rho')
        ndof_per_node : number of DOFs per node (bar element = 1 dof, beam element = 2 dof, 2d beam = 3 dofs)
        ndof_per_node : number of DOFs per node (bar element = 1 dof, beam element = 2 dof, 2d beam = 3 dofs)
    
    Returns:
        K_global : (ndof, ndof) ndarray, the global stiffness matrix
    """

    flag = False
    if rotation_func and ndof_per_node == 1:    # this means we have a 2D truss structure
        ndof_per_node = 2
        flag = True
        

    total_nodes = fem_nodes.shape[0]
    
    total_dofs = ndof_per_node * total_nodes
    K_global = np.zeros((total_dofs, total_dofs))

    E = mat_properties['E']

    # Inspect the signature of the provided local stiffness function
    sig = inspect.signature(k_local_func)
    param_names = list(sig.parameters.keys())

    for e, (i, j) in enumerate(fem_elements):   
        # get start and end point of the fem element
        xi, yi = fem_nodes[i]
        xj, yj = fem_nodes[j]

        # get geometrical properties
        L = np.sqrt((xj - xi)**2 + (yj - yi)**2)
        
        A = element_crossSections[e]['A']
        I = element_crossSections[e]['I']

        # Prepare all possible arguments
        available_args = {
            'E': E,
            'A': A,
            'I': I,
            'L': L
        }

        # Filter only the arguments needed for k_local_func
        args_to_pass = [available_args[name] for name in param_names if name in available_args]


        # Call with only the required arguments
        k_local = k_local_func(*args_to_pass)

        # rotation matrix from local to global - adjust ndof_per_node if flag is up (not the best way to manage the logic, possibily to change)
        if flag:
            ndof_per_node = 1

        if rotation_func == None:
            if ndof_per_node == 1:
                R = np.eye(2)
            else:
                R = np.eye(4)
        else:
            if ndof_per_node == 2:  # need to cut away raws and columns 3 and 6
                R = np.delete(rotation_func(xi,yi,xj,yj), [2, 5], axis=0)  # cut rows
                R = np.delete(R, [2, 5], axis=1)

            elif ndof_per_node == 1:
                R = fem_toolbox.elements.rotation_matrix_bar_2d(xi, yi, xj, yj)

            else:
                R = rotation_func(xi, yi, xj, yj)

        # Transform local stiffness to global coordinates
        k_global = R.T @ k_local @ R

        #import sympy as sy
        #print(f"Global matrix for element {e}: {sy.Matrix(k_global)}")

        # map element DOFs to global DOFs
        if flag:
            ndof_per_node = 2

        dofs = dof_map(ndof_per_node, i, j) 

        # assemble into global matrix
        for a in range(ndof_per_node*2):
            for b in range(ndof_per_node*2):
                K_global[dofs[a], dofs[b]] += k_global[a, b]

    return K_global


def assembleM(m_local_func, rotation_func, fem_nodes, fem_elements, element_crossSections, mat_properties, ndof_per_node):
    """
    Assemble the global mass matrix for 2D frame structure.

    Arguments:
        m_local_func : pointer to the function for local mass matrix evaluation
        rotation_func : pointer to the function for local matrix rotation
        fem_nodes : (N, 2) array of node coordinates
        fem_elements : list of (node_i, node_j)
        element_crossSections : list of dicts with 'A' and 'I'
        mat_properties : dict with 'E' and 'rho' (ρ shared across all elements)
        ndof_per_node : number of DOFs per node (bar element = 1 dof, beam element = 2 dof, 2d beam = 3 dofs)

    Returns:
        M_global : (ndof, ndof) ndarray, the global mass matrix
    """
    flag = False
    if rotation_func and ndof_per_node == 1:    # this means we have a 2D truss structure
        ndof_per_node = 2
        flag = True

    total_nodes = fem_nodes.shape[0]
    total_dofs = ndof_per_node * total_nodes

    M_global = np.zeros((total_dofs, total_dofs))

    # Extract density 
    rho = mat_properties['rho']

    

    for e, (i, j) in enumerate(fem_elements):
        xi, yi = fem_nodes[i]
        xj, yj = fem_nodes[j]

        L = np.sqrt((xj - xi)**2 + (yj - yi)**2)
        A = element_crossSections[e]['A']

        # local mass matrix (consistent mass)
        m_local = m_local_func(rho, A, L)

        # rotation matrix from local to global - adjust ndof_per_node if flag is up (not the best way to manage the logic, possibly to change)
        if flag:
            ndof_per_node = 1

        if rotation_func == None:
            if ndof_per_node == 1:
                R = np.eye(2)
            else:
                R = np.eye(4)
        else:
            if ndof_per_node == 2:  # need to cut away raws and columns 3 and 6
                R = np.delete(rotation_func(xi,yi,xj,yj), [2, 5], axis=0)  # cut rows
                R = np.delete(R, [2, 5], axis=1)

            elif ndof_per_node == 1:
                R = fem_toolbox.elements.rotation_matrix_bar_2d(xi, yi, xj, yj)  # 4x4

            else:
                R = rotation_func(xi, yi, xj, yj)

        # global mass matrix for element
        m_global = R.T @ m_local @ R

        # DOF mapping - correct again for bar elements (to be changed this whole logic)
        if flag:
            ndof_per_node = 2

        dofs = dof_map(ndof_per_node, i, j)

        # assembly
        for a in range(ndof_per_node*2):
            for b in range(ndof_per_node*2):
                M_global[dofs[a], dofs[b]] += m_global[a, b]

    return M_global



def static_analysis(K_global, f_ext, bc_nodes, bc_dofs, bc_values, ndof_per_node):
    """
    Solves the static system K * u = f_ext for the displacement vector u, after applying boundary conditions (BCs).

    Args:
        K_global: The global stiffness matrix.
        f_ext: The external force vector.
        bc_nodes: The nodes where boundary conditions are applied.
        bc_dofs: The DOFs (0=horizontal, 1=vertical, 2=rotation) that are prescribed.
        bc_values: The values of the prescribed boundary conditions.
        ndof_per_node : number of DOFs per node (bar element = 1 dof, beam element = 2 dof, 2d beam = 3 dofs)

    Returns:
        u: The displacement vector, containing displacements for all DOFs.
    """
    
    # 1 create the displacement vector 'u' initialized to zeros
    num_dofs = K_global.shape[0]
    u = np.zeros(num_dofs)

    struct2D = any(i >= 1 for i in bc_dofs)
    flag = True
    # if 2D structure with bar elements (1 dof)
    if ndof_per_node == 1 and struct2D:
        ndof_per_node = 2
        flag = False
        #print("flag false")
    
    # 2 apply boundary conditions (set the prescribed displacements in 'u')
    constrained_dof_ids = [node * ndof_per_node + dof for node, dof in zip(bc_nodes, bc_dofs)]
    #print(f"\nConstrained dofs: {constrained_dof_ids}")

    # corrections if using beam or bar elements
    if ndof_per_node == 2 and flag:
        if any(i == 0 for i in bc_dofs):
            raise IOError("If you use beam elements, you cannot give BC on the x direction. The element has no dof in that direction")
        else:
            # shift indexes
            bc_dofs = np.array(bc_dofs)
            bc_dofs -= 1
            constrained_dof_ids = [node * ndof_per_node + dof for node, dof in zip(bc_nodes, bc_dofs)]
            
    
    if ndof_per_node == 1 and flag:
        if any(i == 1 or i == 2 for i in bc_dofs):
            raise IOError("If you use bar elements, you cannot give BC on the y direction or on the rotation. The element has no dofs in those directions")
        else:
            pass    # no shift necessary here
    
    #print(f"bc_nodes: {bc_nodes}\nbc_dofs: {bc_dofs}\nbc_values: {bc_values}")
    #print(f"Boundary conditions array: {constrained_dof_ids}")

    # Apply BC values to the displacement vector
    for dof_id, value in zip(constrained_dof_ids, bc_values):
        u[dof_id] = value
        #print(dof_id)
    
    # 3. Apply boundary conditions to K_global and f_ext by modifying them
    # -> create a reduced system, eliminating the constrained DOFs
    K_reduced = np.delete(K_global, constrained_dof_ids, axis=0)
    K_reduced = np.delete(K_reduced, constrained_dof_ids, axis=1)
    
    f_ext_reduced = np.delete(f_ext, constrained_dof_ids)

    #import sympy as sy
    #display(sy.Matrix(K_reduced))
    #display(sy.Matrix(f_ext_reduced))
    
    # 4. Solve the reduced system for the unknown displacements
    try:
        u_reduced = np.linalg.solve(K_reduced, f_ext_reduced)
    except np.linalg.LinAlgError as e:
        if "Singular matrix" in str(e):
            raise ValueError(
                "Error: The stiffness matrix is singular. "
                "This can occur if more than one bar element is used per truss in a 2D truss structure. "
                "Please ensure that each truss is modeled with a single bar element."
            ) from e
        else:
            raise  # re-raise other unexpected err
    
    # 5. Insert the solved displacements back into the global displacement vector 'u'
    free_dof_ids = list(set(range(num_dofs)) - set(constrained_dof_ids))
    u[free_dof_ids] = u_reduced
    
    return u


from scipy.linalg import eigh  # Symmetric eigenvalue solver

def modal_analysis(K_global, M_global, bc_nodes, bc_dofs, ndof_per_node, num_modes=5, verbose=True):
    """
    Solves the eigenvalue problem for natural frequencies of a 2D frame.

    Parameters:
        K_global : Global stiffness matrix
        M_global : Global mass matrix
        bc_nodes : List of constrained nodes
        bc_dofs : List of constrained dofs (0: u_x, 1: u_y, 2: theta)
        ndof_per_node : number of DOFs per node (bar element = 1 dof, beam element = 2 dof, 2d beam = 3 dofs)
        num_modes : Number of natural frequencies to compute

    Returns:
        frequencies_hz : Natural frequencies in Hz (array of length num_modes)
        mode_shapes : Mode shapes (each column is a mode)
    """
    total_dofs = K_global.shape[0]
    all_dofs = np.arange(total_dofs)

    flag = True
    struct2D = any(i >= 1 for i in bc_dofs)
    if ndof_per_node == 1 and struct2D:
        ndof_per_node = 2
        flag = False

    # Build constrained DOFs list
    constrained_dofs = [node * ndof_per_node + dof for node, dof in zip(bc_nodes, bc_dofs)]


    # corrections if usingeam or bar elements
    if ndof_per_node == 2 and flag:
        if any(i == 0 for i in bc_dofs):
            raise IOError("If you use beam elements, you cannot give BC on the x direction. The element has no dof in that direction")
        else:
            # shift indexes
            bc_dofs = np.array(bc_dofs)
            bc_dofs -= 1
            constrained_dofs = [node * ndof_per_node + dof for node, dof in zip(bc_nodes, bc_dofs)]


    if ndof_per_node == 1 and flag:
        if any(i == 1 or i == 2 for i in bc_dofs):
            raise IOError("If you use bar elements, you cannot give BC on the y direction or on the rotation. The element has no dofs in those directions")
        else:
            pass    # no shift necessary here

    constrained_dofs = sorted(set(constrained_dofs))

    # Free DOFs
    free_dofs = np.setdiff1d(all_dofs, constrained_dofs)

    # Reduce K and M
    K_ff = K_global[np.ix_(free_dofs, free_dofs)]
    M_ff = M_global[np.ix_(free_dofs, free_dofs)]

    # Solve the generalized eigenvalue problem
    eigvals, eigvecs = eigh(K_ff, M_ff)

    # Remove negative or zero eigenvalues (numerical issues)
    positive = eigvals > 1e-8
    eigvals = eigvals[positive]
    eigvecs = eigvecs[:, positive]

    # Limit to num_modes
    eigvals = eigvals[:num_modes]
    eigvecs = eigvecs[:, :num_modes]

    # Convert to frequencies in Hz
    omegas = np.sqrt(eigvals)
    frequencies_hz = omegas / (2 * np.pi)

    if verbose == True:
        for i, f in enumerate(frequencies_hz):
            print(f"Mode {i+1}: {f:.2f} Hz")


    return frequencies_hz, eigvecs, free_dofs




def build_force_vector(f_nodes, f_dofs, f_values, num_dofs, dofs_per_node, struct2D_trussElements=False):
    """
    Assembles the global external force vector for the structure.

    Args:
        f_nodes: List of node indices where forces are applied.
        f_dofs: List of DOF indices (0: u, 1: v, 2: θ) at each node where force is applied.
        f_values: List of force or moment values to apply.
        num_dofs: Total number of degrees of freedom.
        dofs_per_node: Number of degrees of freedom per node (e.g. 1 for rod, 2 for beam, 3 for 2D frame).

    Returns:
        f_ext: Global force vector (1D numpy array).
    """
    # managing case 1 dof for 2D geometry
    if dofs_per_node == 1 and struct2D_trussElements:
        dofs_per_node = 2

    f_ext = np.zeros(num_dofs)

    for node, dof, value in zip(f_nodes, f_dofs, f_values):
        if dofs_per_node == 2  and not struct2D_trussElements: # if beam element (2 dofs), no axial force can be applied, so shift by 1:
            dof = dof - 1
        global_dof = node * dofs_per_node + dof
        f_ext[global_dof] += value  # add in case multiple forces on same DOF

    return f_ext






