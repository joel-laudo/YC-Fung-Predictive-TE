from abaqus import *
import sys
from odbAccess import *
import numpy as np
from abaqusConstants import *
from visualization import *

# Integration point coordinates at first time step
pts = np.loadtxt('P14_skin_coords_cutout_FINAL.txt')
pts = np.array(pts)
pts_len = len(pts)
pts = pts.astype(np.float32)
increments = np.linspace(11590, 11590, 1)
increments = increments.astype(np.int32)

num_sim = 11
for tmp in range(num_sim):
    tmp = tmp + 170
    if tmp == 6 or tmp == 22 or tmp == 55 or tmp == 71 or tmp == 79:
        tmp = tmp
    else:
        if tmp < 190:
            tmp = tmp
        else:
            tmp = -1
        odbPath = 'P14ExpJobMAP.odb'

        o1 = session.openOdb(name=odbPath)
        odb = session.odbs[odbPath]
        session.viewports['Viewport: 1'].setValues(displayedObject=o1)

        for increment in increments:

            session.viewports['Viewport: 1'].odbDisplay.setFrame(
                step=0, frame=increment)
            session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(
                averageElementOutput=False)
            session.Path(name='Path-1', type=POINT_LIST, expression=(pts))
            pth = session.paths['Path-1']

            session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                variableLabel='COORD', refinement=[COMPONENT, 'COOR1'], outputPosition=NODAL)

            XcoordsWithPath = session.XYDataFromPath(name='XYData-1', path=pth, includeIntersections=False, projectOntoMesh=False,
                                                     pathStyle=PATH_POINTS, numIntervals=10, projectionTolerance=0, shape=UNDEFORMED, labelType=TRUE_DISTANCE)

            Xcoords = [tup[1] for tup in XcoordsWithPath]
            Xcoords = np.array(Xcoords)
            #Xcoords = Xcoords.transpose()

            del session.xyDataObjects['XYData-1']

            session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                variableLabel='COORD', refinement=[COMPONENT, 'COOR2'], outputPosition=NODAL)
            YcoordsWithPath = session.XYDataFromPath(name='XYData-1', path=pth, includeIntersections=False, projectOntoMesh=False,
                                                     pathStyle=PATH_POINTS, numIntervals=10, projectionTolerance=0, shape=UNDEFORMED, labelType=TRUE_DISTANCE)

            Ycoords = [tup[1] for tup in YcoordsWithPath]
            Ycoords = np.array(Ycoords)
            #Ycoords = Ycoords.transpose()

            del session.xyDataObjects['XYData-1']

            session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
                variableLabel='COORD', refinement=[COMPONENT, 'COOR3'], outputPosition=NODAL)
            ZcoordsWithPath = session.XYDataFromPath(name='XYData-1', path=pth, includeIntersections=False, projectOntoMesh=False,
                                                     pathStyle=PATH_POINTS, numIntervals=10, projectionTolerance=0, shape=UNDEFORMED, labelType=TRUE_DISTANCE)

            Zcoords = [tup[1] for tup in ZcoordsWithPath]
            Zcoords = np.array(Zcoords)
            #Zcoords = Zcoords.transpose()

            DeformedCoords = np.array([Xcoords, Ycoords, Zcoords])
            DeformedCoords = DeformedCoords.transpose()
            np.savetxt('P14_Deformed_Coords_Sim_' +
                       str(tmp+1) + '_MAP.txt', DeformedCoords)
        session.odbs[odbPath].close()
