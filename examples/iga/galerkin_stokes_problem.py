r"""
In this file we will prototype a galerkin example (in the simplest way),
that shows how to implement a Stokes problem in the form

.. math::
        -\mu\Delta v + \nabla p = 0

"""
import matplotlib.pyplot as plt

import numpy as np

import splinepy as sp

try:
    import gustaf as gus

    has_gus = True
except ImportError:
    has_gus = False

# Constants for the remainder of this test case
N_REFINE = 25
VISCOSITY = 1e-4
SHOW_MATRIX = True

# Auxiliary function to map positions into the element
def map_positions(positions, x_min, x_max, y_min, y_max):
    mapped_positions = np.empty(positions.shape, dtype=np.float64)
    mapped_positions[:, 0] = (x_max - x_min) * 0.5 * positions[:, 0] + (
        x_max + x_min
    ) * 0.5
    mapped_positions[:, 1] = (y_max - y_min) * 0.5 * positions[:, 0] + (
        y_max + y_min
    ) * 0.5
    return mapped_positions

# Define the geometry (i.e., a unit square for simplicity)
geometry = sp.BSpline(
    degrees=[1, 1],
    control_points=[
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [1.0, 1.0],
    ],
    knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1]],
)

# Define the solution fields
velocity_field = sp.BSpline(
    degrees=geometry.degrees,
    control_points=np.ones((geometry.control_points.shape[0], 2)),
    knot_vectors=geometry.knot_vectors,
)
pressure_field = sp.BSpline(
    degrees=geometry.degrees,
    control_points=np.ones((geometry.control_points.shape[0], 1)),
    knot_vectors=geometry.knot_vectors,
)

# Use refinement
velocity_field.elevate_degrees([0,0, 1,1])
new_knots = np.linspace(1 / N_REFINE, 1, N_REFINE, endpoint=False)
velocity_field.insert_knots(0, new_knots)
velocity_field.insert_knots(1, new_knots)
pressure_field.elevate_degrees([0, 1])
new_knots = np.linspace(1 / N_REFINE, 1, N_REFINE, endpoint=False)
pressure_field.insert_knots(0, new_knots)
pressure_field.insert_knots(1, new_knots)

# Retrieve integration points and weights
max_order = int(np.max([velocity_field.degrees,pressure_field.degrees]))
positions, weights = np.polynomial.legendre.leggauss(deg=max_order)
positions = sp.utils.data.cartesian_product(
    [positions for _ in range(geometry.para_dim)]
)
weights = sp.utils.data.cartesian_product(
    [weights for _ in range(geometry.para_dim)]
)
weights = np.prod(weights, axis=1)

# Get offsets for organizing the degrees of freedom
velocity_offset = velocity_field.control_points.size
vel_dim_offset = velocity_field.control_points.shape[0]

# Determine number of degrees of freedom
velocity_dofs = velocity_offset
pressure_dofs = pressure_field.control_points.size
n_dofs = velocity_dofs + pressure_dofs
print(f'Velocity DoFs: {velocity_dofs}')
print(f'Pressure DoFs: {pressure_dofs}')
print(f'Total DoFs   : {n_dofs}')

# Allocate the resources for solving the linear system
system_matrix = np.zeros((n_dofs, n_dofs))
system_rhs = np.zeros(n_dofs)

# Assemble individual elements
ukv = velocity_field.unique_knots
ukp = pressure_field.unique_knots
velocity_mapper = velocity_field.mapper(reference=geometry)
pressure_mapper = pressure_field.mapper(reference=geometry)

print()
print('Assembling system matrix...', end='', flush=True)

# Element Loop
for i in range(len(ukv[0]) - 1): # iterate over all test velocity DOFs
    for j in range(len(ukv[1]) - 1): # iterate over all trial velocity DOFs
        
        # Update position of current element
        mapped_positions = map_positions(
            positions,
            ukv[0][i],
            ukv[0][i + 1],
            ukv[1][j],
            ukv[1][j + 1],
        )
        # Compute transformation factor
        det_jacs = np.linalg.det(geometry.jacobian(mapped_positions))
        
        # Retrieve gradient and support of the velocity basis functions at the 
        # current position
        vel_bf_grad, vel_bf_sup = velocity_mapper.basis_gradient_and_support(
            mapped_positions
        )
        assert np.all(vel_bf_sup[0, :] == vel_bf_sup) # TODO Was macht das?
        
        # number of velocity basis functions on current element
        num_vel_bf = vel_bf_grad.shape[1]
        # blow up to tensorial shape
        vel_bf_grad_vec = np.einsum(
            "qix,qjy,ik,jk->qkxy",
            vel_bf_grad,
            vel_bf_grad,
            np.eye(num_vel_bf),
            np.eye(num_vel_bf)
        )
        
        ##############################
        # Assemble grad(v) : grad(w) #
        ##############################
        # q : quadrature point | i/j: test/trial function | x,y: dim
        local_matrix = np.einsum(
            "qixy,qjxy,q,q->ij",
            vel_bf_grad_vec,
            vel_bf_grad_vec,
            weights,
            det_jacs,
            optimize=True,
        )

        # construct cartesian structure for mapping local contributions into 
        # global system matrix
        local_to_global = sp.utils.data.cartesian_product(
            (vel_bf_sup[0, :], vel_bf_sup[0, :])
        )
        # write into system matrix
        system_matrix[local_to_global[:, 0], local_to_global[:, 1]] \
            += VISCOSITY*local_matrix.flatten()
        
        ################################
        # Assemble grad(v)^T : grad(w) #
        ################################
        # q : quadrature point | i/j: test/trial function | x,y: dim
        local_matrix = np.einsum(
            "qiyx,qjxy,q,q->ij",
            vel_bf_grad_vec,
            vel_bf_grad_vec,
            weights,
            det_jacs,
            optimize=True,
        )

        # write into system matrix
        system_matrix[local_to_global[:, 0], local_to_global[:, 1]] \
            += VISCOSITY*local_matrix.flatten()
        
    # end iterate over all trial velocity DOFs
    for j in range(len(ukp[1]) - 1): # iterate over all trial pressure DOFs

        # Update position of current element
        mapped_positions = map_positions(
            positions,
            ukv[0][i],
            ukv[0][i + 1],
            ukp[1][j],
            ukp[1][j + 1],
        )
        # Compute transformation factor
        det_jacs = np.linalg.det(geometry.jacobian(mapped_positions))

        # Retrieve gradient and support of the velocity basis functions at the  
        # current position
        vel_bf_grad, vel_bf_sup = velocity_mapper.basis_gradient_and_support(
            mapped_positions
        )
        assert np.all(vel_bf_sup[0, :] == vel_bf_sup) # TODO Was macht das?

        # number of velocity basis functions on current element
        num_vel_bf = vel_bf_grad.shape[1]
        # blow up to tensorial shape
        vel_bf_grad_vec = np.einsum(
            "qix,qjy,ik,jk->qkxy",
            vel_bf_grad,
            vel_bf_grad,
            np.eye(num_vel_bf),
            np.eye(num_vel_bf)
        )
        # Retrieve value and support of the pressure basis functions at the 
        # current position
        pres_bf, pres_bf_sup = pressure_field.basis_and_support(
            mapped_positions
        )
        assert np.all(pres_bf_sup[0, :] == pres_bf_sup) # TODO Was macht das?

        #####################################
        # Assemble p tr(grad(w)) = p div(w) #
        #####################################
        # q : quadrature point | i/j: test/trial function | d: dim
        local_matrix = np.einsum(
            "qj,qidd,q,q->ij",
            pres_bf,
            vel_bf_grad_vec,
            weights,
            det_jacs,
            optimize=True,
        )

        # account for offset by velocity as pressure is the second unkown
        pres_bf_sup += velocity_offset
        # construct cartesian structure for mapping local contributions into 
        # global system matrix
        local_to_global = sp.utils.data.cartesian_product(
            (vel_bf_sup[0, :], pres_bf_sup[0, :])
        )
        # write into system matrix
        system_matrix[local_to_global[:, 0], local_to_global[:, 1]] += local_matrix.flatten()
        
        #####################################
        # Assemble q tr(grad(v)) = q div(u) #
        #####################################
        # this is exactly the transpose of the local matrix assembled above
        local_matrix = np.transpose(local_matrix)

        # write into system matrix
        system_matrix[local_to_global[:, 1], local_to_global[:, 0]] += local_matrix.flatten()
    # end iterate over all trial pressure DOFs
# end iterate over all test velocity DOFs

print('done')

if SHOW_MATRIX:
    plt.spy(system_matrix)
    plt.show()