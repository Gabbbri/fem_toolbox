# Stiffness module
# Definition of local stiffness matrices 

import numpy as np


def k_rod(E, A, L):
    """2-node rod element stiffness matrix"""
    k_local = (E*A/L) * np.array([[1, -1], 
                            [-1, 1]])
    return k_local

def k_beam(E, I, L):
    """2-node beam element stiffness matrix (Euler Bernoulli beam)"""
    k_local = (E*I/L**3) * np.array([[12, 6*L, -12, 6*L],
                             [6*L, 4*L**2, -6*L, 2*L**2],
                             [-12, -6*L, 12, -6*L],
                             [6*L, 2*L**2, -6*L, 4*L**2]]) 
    return k_local

def k_beam2d(E, A, I, L):
    """
    Local stiffness matrix (6x6) for a 2D beam (plane frame) element.
    Includes axial and bending contributions.
    
    DOFs: [u1, v1, θ1, u2, v2, θ2]  (2 nodes × 3 DOFs each)
    
    Parameters:
        E : float  - Young's modulus
        A : float  - Cross-sectional area
        I : float  - Moment of inertia
        L : float  - Length of the element
    
    Returns:
        K_local : (6x6) np.ndarray - Local stiffness matrix
    """
    # Axial stiffness component (rod behavior)
    k_axial = (E * A / L) * np.array([
        [ 1,  0,  0, -1,  0,  0],
        [ 0,  0,  0,  0,  0,  0],
        [ 0,  0,  0,  0,  0,  0],
        [-1,  0,  0,  1,  0,  0],
        [ 0,  0,  0,  0,  0,  0],
        [ 0,  0,  0,  0,  0,  0]
    ])

    # Bending stiffness component (Euler–Bernoulli beam behavior)
    coeff = E * I / L**3
    k_bending = coeff * np.array([
        [ 0,     0,      0,      0,     0,      0     ],
        [ 0,    12,     6*L,     0,   -12,     6*L    ],
        [ 0,   6*L,   4*L**2,    0,  -6*L,   2*L**2   ],
        [ 0,     0,      0,      0,     0,      0     ],
        [ 0,   -12,   -6*L,      0,    12,   -6*L     ],
        [ 0,   6*L,   2*L**2,    0,  -6*L,   4*L**2   ]
    ])

    # Total local stiffness matrix
    k_local = k_axial + k_bending

    return k_local



def m_rod(rho, A, L):
    """Consistent mass matrix for bar element"""
    m_local = rho * A * L / 6 * np.array([
        [2, 1],
        [1, 2]
    ])
    return m_local


def m_beam(rho, A, L):
    """Consistent mass matrix for 2-node Euler-Bernoulli beam element"""
    m_local = rho * A * L / 420 * np.array([
        [156, 22*L, 54, -13*L],
        [22*L, 4*L**2, 13*L, -3*L**2],
        [54, 13*L, 156, -22*L],
        [-13*L, -3*L**2, -22*L, 4*L**2]
    ])
    return m_local


def m_beam2d(rho, A, L):
    """
    Local consistent mass matrix (6x6) for a 2D beam (plane frame) element.
    Includes both axial and bending contributions.
    
    DOFs: [u1, v1, θ1, u2, v2, θ2]

    Parameters:
        rho : float  - Material density (kg/mm³)
        A   : float  - Cross-sectional area (mm²)
        L   : float  - Length of the element (mm)

    Returns:
        M_local : (6x6) np.ndarray - Local mass matrix
    """

    m = rho * A * L  # total mass of element

    # Axial (rod) mass terms — 2x2 consistent mass matrix for u DOFs
    m_axial = (m / 6) * np.array([
        [2, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0],
        [1, 0, 0, 2, 0, 0],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0]
    ])

    # Bending mass matrix from Euler–Bernoulli theory (for v and theta)
    coeff = rho * A * L / 420
    m_bending = coeff * np.array([
        [0,    0,     0,    0,    0,     0],
        [0,   156,   22*L,   0,   54,   -13*L],
        [0,  22*L,  4*L**2,  0, 13*L,  -3*L**2],
        [0,    0,     0,    0,    0,     0],
        [0,    54,   13*L,   0,  156,  -22*L],
        [0, -13*L, -3*L**2,  0, -22*L, 4*L**2]
    ])

    # Total local mass matrix
    m_local = m_axial + m_bending

    return m_local


def rotation_2d(x1, y1, x2, y2):
    """
    Returns the 6x6 rotation matrix for a 2D beam element
    between points (x1, y1) and (x2, y2).

    Parameters:
        x1, y1: coordinates of the first node
        x2, y2: coordinates of the second node

    Returns:
        T: (6x6) np.ndarray - Rotation matrix
    """
    # Direction cosines
    dx = x2 - x1
    dy = y2 - y1
    L = np.sqrt(dx**2 + dy**2)

    if L == 0:
        raise ValueError("Element length is zero. Check node coordinates.")

    c = dx / L
    s = dy / L

    # Rotation matrix R 
    R = np.zeros((6, 6))    # every node has 3 dofs

    # Fill sub-blocks
    R_local = np.array([
        [c, s, 0],
        [-s, c, 0],
        [0, 0, 1]
    ])

    # Upper-left 3x3
    R[0:3, 0:3] = R_local
    # Lower-right 3x3
    R[3:6, 3:6] = R_local

    return R


def rotation_matrix_bar_2d(xi, yi, xj, yj):
        L = np.sqrt((xj - xi)**2 + (yj - yi)**2)
        c = (xj - xi) / L
        s = (yj - yi) / L
        R = np.array([
        [c, s, 0, 0],
        [0, 0, c, s]
        ])
        return R
