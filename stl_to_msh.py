
# coding: utf-8

# # .stl to .msh Python script.
# 
# The goal of this script is to automatize the meshing of a .stl dataset.
# The meshing will be done using *gmsh* library.
# follow this link to get info concerning *gmsh*:https://gmsh.info/doc/texinfo/gmsh.html#index-Mesh_002eQualityType

'''
Import the libraries
'''

import datetime
import pathlib
import time
import multiprocessing as mp
import os
import gmsh


import argparse

parser = argparse.ArgumentParser(description='Process stl files to msh files.')
parser.add_argument("-p" , '--path_stl', metavar='P', type=str, nargs='+',
                    help='Where are stored the stl files')
args = parser.parse_args()
'''
Mesh Parameters
'''


# we have a class where the parameters are

class MeshParameters:
    def __init__(option, Algorithm3D, MinLength, MaxLength, OptimizationNet,
                 Quality, Version, Iteration, AngleTol):
        # kind of mesh, 1 for Delaunay meshes
        option.Algorithm3D = gmsh.option.setNumber("Mesh.Algorithm3D", Algorithm3D)
        option.MinLength = gmsh.option.setNumber("Mesh.CharacteristicLengthMin", MinLength)
        option.MaxLength = gmsh.option.setNumber("Mesh.CharacteristicLengthMax", MaxLength)
        # Optimize the mesh using Netgen to improve the quality of tetrahedral elements
        option.OptimizationNet = gmsh.option.setNumber("Mesh.OptimizeNetgen", OptimizationNet)
        # type of quality 2: gamma~vol/sum_face/max_edge
        option.Quality = gmsh.option.setNumber("Mesh.QualityType", Quality)
        # version of mesh file, use 2.2
        option.Version = gmsh.option.setNumber("Mesh.MshFileVersion", Version)
        # Maximum number of point insertion iterations in 3D Delaunay
        option.Iteration = gmsh.option.setNumber("Mesh.MaxIterDelaunay3D", Iteration)
        # Consider connected facets as overlapping
        option.AngleTol = gmsh.option.setNumber("Mesh.AngleToleranceFacetOverlap", AngleTol)


'''
Report function
'''


def report(file, problem, report_file):
    # we get the problems inside a list
    t = datetime.datetime.now()
    pb_list = [t, ' The file ', file, " has the problem: ", problem, '\n']
    with open(report_file, "a") as fileout:
        for pb in pb_list:
            # we write the problem on the text file
            fileout.writelines(format(pb))


'''
Mesh function
'''

the_lock = mp.Lock()



def mesher(the_stl ):
    # --PATH-- #
    folder = args.path_stl[0]
    if folder[-1] != "/":
        folder = folder + "/"
    # Put your own path:
    folder_stl = folder
    folder_msh = folder
    folder_report = folder + "report.txt"

    # --MANAGEMENT-- #

    # we # print the CPU working on the stl
    # print('The', the_stl, 'is on', mp.current_process())
    # Manage the name
    x = the_stl.split('.')
    file_name = x[0]

    # --MESHING-- #

    # print("Starting job at", datetime.datetime.now())
    gmsh.initialize()
    gmsh.clear()

    # we fill the class with parameters
    MeshParameters(1,  # Algorithm
                   0.01,  # MinLength
                   15.00,  # MaxLength
                   0,  # Optimization Net
                   2,  # Quality
                   2.2,  # Version
                   2155000,  # iteration
                   0.05)  # angle tolerance


    try:
        # print("Processing", the_stl)
        gmsh.merge(folder_stl + the_stl)  # we capture the stl file

        n = gmsh.model.getDimension()
        s = gmsh.model.getEntities(n)
        surf = gmsh.model.geo.addSurfaceLoop([s[i][1] for i in range(len(s))])

        gmsh.model.geo.addVolume([surf])  # Create one volume
        # print("Volume added")
        gmsh.model.geo.synchronize()
        # dimen = gmsh.model.getDimension()
        # we generate the mesh here
        gmsh.model.mesh.generate(3)
        # when  the generation is done we write the .msh file
        gmsh.write(folder_msh + file_name + '.msh')  # we create the .msh file

        gmsh.finalize()

        # print("Meshing is done")
        # print("")


    except Exception as exception:
        # print("  _/!\_File", the_stl, 'has a problem. Problem:', exception)
        report(the_stl, exception, folder_report)

        gmsh.finalize()

        # we manage the exception here
        # if str(exception)== 'Invalid boundary mesh (overlapping facets) on surface 1 surface 1':
        # gmsh.model.geo.remove([(1,1)]) #we solve the problem
        # mesh_here(folder_msh,fileName)#we generate again
        # # print("")

        # if str(exception) == 'PLC Error:  A segment and a facet intersect at point':
        # gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        # mesh_here(folder_msh,fileName)#we generate again
        # print("")


def mesher2(the_stl:str, folder_stl:str, folder_msh:str, folder_report:str):

    # --MANAGEMENT-- #
   
    # Manage the name
    x = the_stl.split('.')
    file_name = x[0]

    # --MESHING-- #

    # print("Starting job at", datetime.datetime.now())
    # print(0)
    gmsh.initialize()
    # print(1)
    gmsh.clear()
    # print(2)

    # we fill the class with parameters
    MeshParameters(1,  # Algorithm
                   0.01,  # MinLength
                   15.00,  # MaxLength
                   0,  # Optimization Net
                   2,  # Quality
                   2.2,  # Version
                   2155000,  # iteration
                   0.05)  # angle tolerance


    try:
        # print("Processing", the_stl)
        gmsh.merge(folder_stl + the_stl)  # we capture the stl file

        n = gmsh.model.getDimension()
        s = gmsh.model.getEntities(n)
        surf = gmsh.model.geo.addSurfaceLoop([s[i][1] for i in range(len(s))])

        gmsh.model.geo.addVolume([surf])  # Create one volume
        # print("Volume added")
        gmsh.model.geo.synchronize()
        # dimen = gmsh.model.getDimension()
        # we generate the mesh here
        gmsh.model.mesh.generate(3)
        # when  the generation is done we write the .msh file
        gmsh.write(folder_msh + file_name + '.msh')  # we create the .msh file

        gmsh.finalize()

        # print("Meshing is done")
        # print("")


    except Exception as exception:
        # print("  _/!\_File", the_stl, 'has a problem. Problem:', exception)
        report(the_stl, exception, folder_report)

        gmsh.finalize()

        # we manage the exception here
        # if str(exception)== 'Invalid boundary mesh (overlapping facets) on surface 1 surface 1':
        # gmsh.model.geo.remove([(1,1)]) #we solve the problem
        # mesh_here(folder_msh,fileName)#we generate again
        # print("")

        # if str(exception) == 'PLC Error:  A segment and a facet intersect at point':
        # gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        # mesh_here(folder_msh,fileName)#we generate again
        # # print("")

'''
Check function
'''


def checker(stl_folder):
    # print(stl_folder)
    folder = os.listdir(stl_folder)
    with os.scandir(stl_folder) as ENTRIES:
        c = 0
        for file in ENTRIES:
            # Inside the folder with .stl, we catch each .stl file to mesh them.
            if file.name.split('.')[1] == 'stl':
                c = c + 1
    # print("There are", c, ".stl files inside the dataset.")
    # print()



'''
MAIN
'''


if __name__ == '__main__':
    # --PATH-- #

    path_stl = args.path_stl[0]

    if path_stl[-1] != "/":     
        path_stl = path_stl + "/"

    # --CHECK-- #

    checker(path_stl)

    # --MULTIPROCESSING PARAMETERS-- #

    first_time = datetime.datetime.now()
    # print('Number of CPUs available:', mp.cpu_count())

    # Put inside a list the stl that have to be meshed
    big_list = [] # all the stl here
    small_list = [] # just a small number of stl to do not have memory leak
    with os.scandir(path_stl) as entries:

        for f in entries:
            if f.name.split('.')[1] == 'stl':
                big_list.append(str(f.name))

        big_size = len(big_list)
        # print(big_size)

        while len(big_list) != 0:
            element = 8
            # here the collect some element from the big list
            small_list = big_list[:element]
            del big_list[:element]
            # print("")
            # print('----THE SMALL given to the mp :', small_list)
            # print("")

            # --MESH WITH MULTIPROCESSING-- #
            pool = mp.Pool(processes=8)  # define number of cpus
            pool.map(mesher, small_list)

            #for i in range(len(small_list)):
            #    os.remove(path_stl+'/'+small_list[i])

            # --MANAGE LIST FOR MP-- #

            del small_list[:element]
            pool.close()
            pool.join()
            pool.terminate()
            # print('When terminate, big:', len(big_list), ' small', len(small_list))
            time.sleep(3)

    later_time = datetime.datetime.now()
    # print("")
    # print("---- Total time: %s " % (later_time - first_time))
    # print("-- Dataset meshing is done :) --")
