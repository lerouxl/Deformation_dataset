import logging
from tqdm import tqdm
from pathlib import Path
from InherentStrain import InherentStrain
from cubes_generator import generate_cubes
from check_edges_number import remove_not_exact
import os
import pyvista as pv
import numpy as np
import trimesh 

logger_path = "Logs.log"
logging.basicConfig(format='%(asctime)s - %(message)s', filename=logger_path, level=logging.DEBUG)

def generate_dataset(dataset_path: Path, final_path: Path, number_sample: int = 100,
                    max_edge_size: float = 2.0, noise: float = 0.001, edges_target: int = 3000,
                    number_of_vert: int = 2000, gmsh_path: str = r"gmsh-3.0.6-Linux64/bin/gmsh" ):
    """
    Generate the dataset for AM simulation.
    Args:
        - dataset_path: Where are save the temp files
        - final_path: Where the dataset is saved
        - number_sample: how many stl files to generate
        - max_edge_size: max size of edges, BEFORE DECIMATION
        - noise: size of the noise added to the cubes
        - edges_target: number of edges to have per file
        - number_of_vert defimation target
        - gmsh_path: Where is gmsh 3.0.6

    """
    # Create the dataset folder if it do not exist
    dataset_path.mkdir(parents=True, exist_ok=True)
    final_path.mkdir(parents=True, exist_ok=True)

    # Generate the cubes obj files
    generate_cubes(path = dataset_path, number_sample = number_sample, noise=noise, number_of_vert= number_of_vert, max_edge_size = max_edge_size)

    # Remove files with the wrong number of edges:
    remove_not_exact(folder_path= dataset_path, edges_target=edges_target )
    # Generate the msh file of each obj
    os.system(f"python3 stl_to_msh.py --path_stl {dataset_path}")

    # Generate the simulation of each msh
    msh_files = dataset_path.glob("*.msh")
    msh_files = list(msh_files)
    logging.info(f"Start computing the deformation, {len(msh_files)} files found")
    pbar = tqdm(msh_files)
    for msh in pbar:
        pbar.set_description("Computing deformation files")
        InherentStrain(input= msh, output_folder= dataset_path, k = 1.0e8, z_clamping_tolerance = 0.1)

    # Transfert the simulation results and mesh to the final_path folder
    # Copy the logs as there is the generation option
    src = Path(logger_path)
    dest = final_path / Path(logger_path).name
    dest.write_bytes(src.read_bytes())

    vtu_files = dataset_path.glob("*_Displacement000000.vtu")
    vtu_files = list(vtu_files)

    logging.info(f"Start exporting the results to the final folder, {len(vtu_files)} files found")
    pbar = tqdm(vtu_files)
    for vtu in pbar:
        pbar.set_description("Creating the final dataset")

        unstructured_mesh = pv.read(vtu)
        mesh = unstructured_mesh.extract_surface()

        # Save the mesh
        mesh.save(final_path / Path(vtu.name.split("._")[0]).with_suffix(".vtk"))
        mesh.save(final_path / Path(vtu.name.split("._")[0]).with_suffix(".stl"), binary=False)
    
        # Save the displacement label
        np.savetxt( final_path / Path(vtu.name.split("._")[0]).with_suffix(".txt") , mesh.point_data["displacement"])

        # Load and the mesh in obj
        mesh = trimesh.load(final_path / Path(vtu.name.split("._")[0]).with_suffix(".stl"), process= True)
        mesh.export(final_path / Path(vtu.name.split("._")[0]).with_suffix(".obj"))



for phase in tqdm(["train", "test", "validation"]):

    # where will be stored the dataset
    dataset_path = Path("cubes_process") / phase # Where all the temps files are saved
    final_path = Path("cubes") / phase # Where the AI will take his inputs

    if phase == "train":
        number_sample = 6000
    if phase == "test":
        number_sample = 2000
    elif phase == "validation":
        number_sample = 500

    max_edge_size = 1.4
    noise = 0.005
    number_of_vert = 2000
    edges_target = 3000
    gmsh_path = r"gmsh-3.0.6-Linux64/bin/gmsh"

    logging.info(f"Start generating dataset")
    logging.info(f"Phase {phase} parameters: \n\tnumber_sample:{number_sample} \n\tnoise:{noise} \n\tnumber_of_vert:{number_of_vert} \n\tedges_target:{edges_target} \n\tmax_edge_size:{max_edge_size}")
    generate_dataset(dataset_path, final_path, number_sample,max_edge_size, noise,edges_target, number_of_vert,gmsh_path)

