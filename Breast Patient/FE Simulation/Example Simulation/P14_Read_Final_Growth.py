from abaqus import *
import sys
from odbAccess import *
import numpy as np
from abaqusConstants import *
from visualization import *

# Integration point coordinates at first time step
pts = np.loadtxt('P14_int_pts_coords.txt')
pts = np.array(pts)
pts_len = len(pts)
pts = pts.astype(np.float32)
increments = np.linspace(11590, 11590, 1)
increments = increments.astype(np.int32)

num_sim = 1
for tmp in range(num_sim):
    tmp = tmp + 161
    if tmp == 6 or tmp == 22 or tmp == 55 or tmp == 71 or tmp == 79:
        tmp = tmp
    else:
        if tmp < 190:
            tmp = tmp
        else:
            tmp = -1
        odbPath = 'P14ExpJobMAP.odb'
        th_data = np.zeros([len(increments), pts_len])
        thg_data = np.zeros([len(increments), pts_len])
        row = 0

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
                variableLabel='SDV1', outputPosition=INTEGRATION_POINT)

            thetag = session.XYDataFromPath(name='XYData-1', path=pth, includeIntersections=False, projectOntoMesh=False,
                                            pathStyle=PATH_POINTS, numIntervals=10, projectionTolerance=0, shape=UNDEFORMED, labelType=TRUE_DISTANCE)

            thetagRawIntPtValues = [tup[1] for tup in thetag]

            thetagElemValues = []
            counter = 0
            tempSum = 0
            # iterate through all integration point values which come out in sets of 8 in the order that
            # elements are defined
            for i in range(len(thetagRawIntPtValues)):
                tempSum += thetagRawIntPtValues[i]
                counter += 1
                if counter == 8:
                    counter = 0
                    thetagElemValues.append(tempSum/8)
                    tempSum = 0

            # Save raw growth at all integration points
            np.savetxt('P14_thg_raw_intPt_data_Sim_' +
                       str(tmp+1) + '_MAP.txt', thetag, fmt='%.6f')

            np.savetxt('P14_thg_elem_data_Sim_' + str(tmp+1) +
                       '_MAP.txt', thetagElemValues, fmt='%.6f')

            del session.xyDataObjects['XYData-1']


        session.odbs[odbPath].close()
