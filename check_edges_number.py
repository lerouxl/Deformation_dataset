from cmath import log
from pathlib import Path
import trimesh
import logging

def remove_not_exact(folder_path: Path, edges_target:int = 3000):
    """Remove the mesh with the wrong number of edges.
    Args:
        - folder_path Path: folder to check for stl files
        - edges_target: int number of edges to have. If the mesh do not have the good number of edges, the mesh will be removed
    """
    for file in Path(folder_path).glob("*.stl"):
        mesh = trimesh.load(file)
        edges_number = mesh.edges_unique.shape[0]
        removed_count = 0

        if edges_number != edges_target:
            # If the number is wrong, delect the obj and the stl
            file.unlink()
            file.with_suffix(".obj").unlink
            removed_count +=1
    
    logging.info(f"{removed_count} files where removed due to incorrect edges numbers (target {edges_target})")