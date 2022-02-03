from dolfin import *
import numpy as np
import os
from msh2xdmf import import_mesh_from_xdmf
from pathlib import Path
from scipy.sparse import csc_matrix # pip install scikit-sparse
from scipy.sparse.linalg import spsolve
import sksparse.cholmod
import scipy
import meshio
import logging

#https://github.com/floiseau/msh2xdmf

def InherentStrain(input: Path, output_folder: Path, k: int = 1.0e8, z_clamping_tolerance: float = 0.1):
    """
    Function deforming a part using the inherent strain methods.
    Args:
        input: pathlib.Path, where the input file is
        output_folder: pathlib.Path, where to save the output files
        k = 1.0e8 int, the support rigidity
        z_clamping_tolerance=0.1 float, all vertices under this z will be considered as supports.
    """

    # Add a variable for the input file
    # Add variable for the output folder
    k=5.0e0 # clamping at z< z_clamping_tolerance, support rigidity
    z_clamping_tolerance = 0.1 # ALL nodes with a z lower than his values will be considered as supports
    tol = 1.0e-10
    E = 1.0
    nu = 0.3
    mu = E/(2.0*(1.0 + nu))
    lmbda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu))
    name = Path(input.stem)

    #https://fenicsproject.discourse.group/t/transitioning-from-mesh-xml-to-mesh-xdmf-from-dolfin-convert-to-meshio/412/2
    #https://fenicsproject.discourse.group/t/pygmsh-tutorial/2506/4
    msh = meshio.read(input)
    logging.info(f"Start simulating {input}")
    for cell in msh.cells:  
        if cell.type == "triangle":
            triangle_cells = cell.data
        elif  cell.type == "tetra":
            tetra_cells = cell.data

    tetra_mesh = meshio.Mesh(points=msh.points, cells={"tetra": tetra_cells})
    xdmf_path =  output_folder / name.with_suffix('.xdmf')
    meshio.write(xdmf_path, tetra_mesh)

    mesh_file = XDMFFile(str(xdmf_path))
    mesh = Mesh()
    mesh_file.read(mesh);

    V = VectorFunctionSpace(mesh, 'P', 2)
    #W = FunctionSpace(mesh, 'P', 1)
    T = TensorFunctionSpace(mesh, 'DG', 0)

    # Define boundary condition

    class Boundary_Bottom_func(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and x[2] < z_clamping_tolerance

    class BoundarySides_x_func(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[0],-31,0.5)

    Boundary_Bottom = Boundary_Bottom_func()
 
    BoundaryMarker = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    BoundaryMarker.set_all(0)
    Boundary_Bottom.mark(BoundaryMarker, 1)

    FileBoundaryMarker = File(str(output_folder / f"{name}._BoundaryMarker.pvd"))
    FileBoundaryMarker << BoundaryMarker

    def sigma(v):
        return lmbda*tr(sym(grad(v)))*Identity(len(v)) + 2.0*mu*sym(grad(v))

    class InherentStrainFunc(UserExpression):
        def eval(self, value, x):
            value[0] = -1.0
            value[1] = 0.
            value[2] = 0.
            value[3] = 0.
            value[4] = -1.0
            value[5] = 0.
            value[6] = 0.
            value[7] = 0.
            value[8] = -0.5


    InherentStrain = InherentStrainFunc(element=T.ufl_element())
    InherentStrainInterp = Function(T)


    InherentStrainInterp.interpolate(InherentStrain)
    InherentStrain_file = File( str(output_folder / f"{name}._InherentStrain.pvd"))
    InherentStrain_file << InherentStrainInterp

    u = TrialFunction(V)
    v = TestFunction(V)

    ds = Measure('ds', domain=mesh, subdomain_data=BoundaryMarker)

    a = inner(sigma(u), grad(v))*dx #+ 1.e-5 * inner(u,v) *dx
    a += k*inner(u, v)*ds(1)

    l = inner(InherentStrain,sym(grad(v))) * dx

    u_sol = Function(V)

    solve(a == l, u_sol)

    file_displacement = File( str(output_folder / f"{name}._Displacement.pvd"))
    
    u_sol.rename('displacement','displacement')
    file_displacement << u_sol
    logging.info(f"End simulating {input}")

if __name__ == "__main__":
    input = Path("cube.msh")
    output = Path("results/")

    InherentStrain(input=input, output_folder=output , k=1e8, z_clamping_tolerance=0.1)