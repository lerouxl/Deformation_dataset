import trimesh
import pyvista as pv
import numpy as np
from pathlib import Path
from tqdm import tqdm
from check_edges_number import remove_not_exact

def extraction_correction(path: Path):
    """
    The previous extraction from the vtk files is wrong, this script have to been launched at the end of the dataset generation to correct the data.

    Arg:
        path: pathlib.Path where are the vtk files stored
    """

    path = Path(path)

    vtk_files = list(path.glob("*.vtk"))
    
    banned_name = ["label_extract", "results"] # If the files name contain those strings, it will be skiped

    for vtk_file in tqdm(vtk_files, total=len(vtk_files)):

        if vtk_file.stem in banned_name:
            continue
        vtk = pv.read(str(vtk_file))
        mesh = trimesh.base.Trimesh(vertices = vtk.points , faces = vtk.faces.reshape(-1, 4)[:,1:], process=False) 
        trimesh.repair.fix_normals(mesh)

        mesh.export( vtk_file.with_suffix(".obj"))
        mesh.export( vtk_file.with_suffix(".stl"))
        #trimesh.repair.fix_normals(mesh)
        #radii = np.linalg.norm( vtk.get_array("displacement"), axis=1) 
        #mesh.visual.vertex_colors = trimesh.visual.interpolate(radii, color_map='viridis')
        #mesh.show()

if __name__ == "__main__":
    for phase in ["validation" , "train", "test"]:
        target =  Path("cubes") / phase
        extraction_correction(target)
        remove_not_exact(folder_path= target, edges_target=3000 )