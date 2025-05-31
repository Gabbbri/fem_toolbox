import numpy as np


def read_structure(file_path):
    """
    Reads a frame definition file where each line is:
        x1 y1 x2 y2 A I  (coordinates of beam endpoints and cross-section properties)
    The first line specifies the material properties: E  rho 0 0 0 0

    Returns:
        node_coords_beam : (N, 2) ndarray - unique node coordinates
        connectivity_beam : (E, 2) ndarray - element connectivity as node indices
        beam_crossSections : list fo dicts per beam with A and I
    """
    raw_data = np.loadtxt(file_path, comments='#')
    
    # simple check on 2nd line
    if raw_data.shape[1] != 6:
        raise ValueError("Each line must contain 4 numbers: x1 y1 x2 y2 A I")

    # material properties
    E = raw_data[0,0]
    rho = raw_data[0,1]
    mat_properties = {"E":E, "rho":rho}

    # extracting coordinates
    coords = raw_data[1:, :4]
    cross_sections_params = raw_data[1:, 4:]

    # Collect all endpoints in one vertical vector (multiple are surely present!)
    all_points = np.vstack([coords[:, :2], coords[:, 2:]])

    # Round to avoid floating-point duplicates (if you generate points with software)
    all_points = np.round(all_points, decimals=5)

    # Cut away the repeated points
    node_coords_beam, _ = np.unique(all_points, axis=0, return_inverse=True)

    # Build connectivity_beam (using the INDEXES of the points)
    num_beams = coords.shape[0]
    connectivity_beam = np.zeros((num_beams, 2), dtype=int)
    
    # array for collecting info on cross section for each beam
    beam_crossSections = []

    for i in range(num_beams):
        p1 = coords[i, :2]
        p2 = coords[i, 2:]
        
        # finding the index related with the specific point
        idx1 = np.where((node_coords_beam == np.round(p1, 5)).all(axis=1))[0][0]
        idx2 = np.where((node_coords_beam == np.round(p2, 5)).all(axis=1))[0][0]

        connectivity_beam[i] = [idx1, idx2]
        beam_crossSections.append({
            "A": cross_sections_params[i, 0],
            "I": cross_sections_params[i, 1]
        })


    return node_coords_beam, connectivity_beam, beam_crossSections, mat_properties


