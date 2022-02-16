from pymeshlab import Mesh as pmMesh
from pymeshlab import MeshSet as pmMeshSet
import trimesh
from pathlib import Path
from tqdm import tqdm
import random
import logging
import numpy as np
from tqdm import tqdm 
import random
import pyvista as pv
import logging
import os

def reduce(mesh:trimesh.base.Trimesh, max_faces_number: int, show_evolution: bool= False):
    """
    Take a trimesh mesh and reduce his number of faces to reach max_faces_number. 
    """
    m = pmMesh(mesh.vertices, mesh.faces)
    del mesh
    ms = pmMeshSet()
    ms.add_mesh(m,"part")

    short_history = [] # If the number of faces to not decrease, we should stop the reduction as there is a probleme
    # If there is too mutch faces:
    if ms.current_mesh().face_number() > max_faces_number:
        # We will reduce the number of face iteratively using a linear equation ax+b
        # with x the step number (we would like to to it in at least 10 step)

        # Number of step to smothely reduce the number of vertice
        number_of_step = (ms.current_mesh().face_number() - max_faces_number) //100
        if number_of_step < 10:
            number_of_step = 10
        elif number_of_step > 10000:
            number_of_step = 10000
        logging.debug(f"\t will be reduced in {number_of_step} steps")

        b = ms.current_mesh().face_number()
        a = (max_faces_number - b) / number_of_step
        x = 0
        count = 0

        # Simplify the mesh. Only first simplification will be agressive
        while (ms.current_mesh().face_number() > max_faces_number):
            if x < number_of_step:
                numFaces = int(a * x + b)
            else:
                numFaces = max_faces_number
                # Avoid error by making sure that the last step had the exact number of faces
            if numFaces < max_faces_number:
                numFaces = max_faces_number

            ms.apply_filter('simplification_quadric_edge_collapse_decimation', targetfacenum=numFaces)
            x += 1
            if x % 100 == 0:
                #print(f"Target {max_faces_number} actual {ms.current_mesh().face_number()} (actual target {numFaces})")

                # If the number of faces to not decrease, we should stop the reduction as there is a probleme
                short_history.append(ms.current_mesh().face_number())
                if len(short_history) > 3:
                    if short_history[-1] == short_history[-2] == short_history[-3]:
                        # if there is no evolution stop
                        raise Exception("The mesh is not decimable with those parameters")

    else:
        msg = " do not have enough faces, skipped"
        raise Exception(msg)

    m = ms.current_mesh()
    vertex_array_matrix = m.vertex_matrix()
    face_array_matrix = m.face_matrix()

    mesh = trimesh.Trimesh(vertices=vertex_array_matrix, faces=face_array_matrix, process=False)
    return mesh


def generate_cubes(path: Path, number_sample: int = 1, noise: float = 0.001, name_prefix: str = "cube_", 
                    height_max: float = 30, height_min: float = 10 ,  max_edge_size: float = 1.0 ,number_of_vert: int = 1500):
    """
    Generate number_samble cubes with random dimension.
    Args:
        path: pathlib.Path: where to save the result
        number_sample: int, how many cubes to generate
        noise: float: size of the noise to add at each cube. This cubes his needed for the faces reducion step.
        name_prefix: str, prefix added to the name of the saved cubes.
        height_max: float = 30 Maximum size of the cubes (not just height, deep and width of the cubes)
        height_min: float = 10 Minimum size of the cubes (not just height, deep and width of the cubes)
        max_edge_size: float = 1.0 Maximum size of the edges BEFORE the decimation
        number_of_vert: int = 1500 Number of faces AFTER the decimation (will change the max edge size of the mesh)
    """
    for i in tqdm(range(number_sample), total=number_sample):
        try:
            logging.debug(f"cubes N {i}: starting dimension generation")
            z_lenght = random.uniform(height_min, height_max)
            x_lenght = random.uniform(height_min, height_max)
            y_lenght = random.uniform(height_min, height_max)
            box = trimesh.creation.box(extents=(x_lenght, y_lenght, z_lenght))

            vert , faces = trimesh.remesh.subdivide_to_size(box.vertices, box.faces, max_edge=max_edge_size)
            logging.debug(f"cubes N {i}: subdivision finished")

            box = trimesh.base.Trimesh(vertices =vert  , faces = faces)
            box = trimesh.permutate.noise(box, noise)
            logging.debug(f"cubes N {i}: noised added")
            box = reduce(box, number_of_vert)
            logging.debug(f"cubes N {i}: reduce finished")

            # Put z_min = z_mov = 0
            z_min = - min(box.vertices[:,2])
            z_mov = trimesh.transformations.translation_matrix(direction=[0,0,z_min])
            box.apply_transform(z_mov)
            export_path = path / f"{i}.obj"     
            box.export(export_path)
            box.export(export_path.with_suffix(".stl"))
            # remove the first line of the obj, with is a comment and made GMSH crash
            #os.system(f"sed -i '1d' {export_path}")
            logging.debug(f"cubes N {i}: Finished")
        except Exception as e:
            logging.critical(e, exc_info=True)

def generate_polygon(path: Path, number_sample: int = 1, noise: float = 0.001, name_prefix: str = "cube_", 
                    height_max: float = 30, height_min: float = 10 ,  max_edge_size: float = 1.0 ,number_of_vert: int = 1500):
    """
    Generate number_samble polygon with random dimension and number of edges.
    Args:
        path: pathlib.Path: where to save the result
        number_sample: int, how many cubes to generate
        noise: float: size of the noise to add at each cube. This cubes his needed for the faces reducion step.
        name_prefix: str, prefix added to the name of the saved cubes.
        height_max: float = 30 Maximum size of the cubes (not just height, deep and width of the cubes)
        height_min: float = 10 Minimum size of the cubes (not just height, deep and width of the cubes)
        max_edge_size: float = 1.0 Maximum size of the edges BEFORE the decimation
        number_of_vert: int = 1500 Number of faces AFTER the decimation (will change the max edge size of the mesh)
    """
    for i in tqdm(range(number_sample), total=number_sample):
        try:
            logging.debug(f"cubes N {i}: starting dimension generation")
            height = random.uniform(height_min, height_max)
            segments = random.randint(4,10)
            radius = random.uniform(height_min, height_max)
            polygon = trimesh.path.polygons.random_polygon(segments= segments, radius = radius)
            extrusion = trimesh.primitives.Extrusion(polygon=polygon, height = height)
            vert , faces = trimesh.remesh.subdivide_to_size(extrusion.vertices, extrusion.faces, max_edge=max_edge_size)

            logging.debug(f"cubes N {i}: subdivision finished")

            box = trimesh.base.Trimesh(vertices =vert  , faces = faces)
            box = trimesh.permutate.noise(box, noise)
            logging.debug(f"cubes N {i}: noised added")
            box = reduce(box, number_of_vert)
            logging.debug(f"cubes N {i}: reduce finished")

            # Put z_min = z_mov = 0
            z_min = - min(box.vertices[:,2])
            z_mov = trimesh.transformations.translation_matrix(direction=[0,0,z_min])
            box.apply_transform(z_mov)
            export_path = path / f"{i}.obj"     
            box.export(export_path)
            box.export(export_path.with_suffix(".stl"))
            # remove the first line of the obj, with is a comment and made GMSH crash
            #os.system(f"sed -i '1d' {export_path}")
            logging.debug(f"cubes N {i}: Finished")
        except Exception as e:
            logging.critical(e, exc_info=True)
