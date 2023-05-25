# -*- coding: utf-8 -*-
"""
Created on Wed Aug 02 19:16:16 2017

@author: G. Zavodszky
@note: Requires vtk
"""

import random
import math
import time
import sys
import vtk  # The python wrapper of VTK

            
class SpatialNode(object):
    def __init__(self, ID, time):
        self.ID = ID
        self.t = time
        self.future_nodes = []
        self.past_nodes = []
        self.right_neighbour = None 
        self.boundary = 0  

class Universe(object):
    def __init__ (self):
        self.spatialNodes = {}
        self.numTimeSlices = 0
        self.nextNodeID = -1
        self.scaledCosmologicalConstant = 0.693
        self.boltzmannWeight_neg = math.exp(-self.scaledCosmologicalConstant)
        self.boltzmannWeight_pos = math.exp(self.scaledCosmologicalConstant)
        self.moveCounter = 0
        self.invMoveCounter = 0

        # Visualization stuff
        self.render = None
        self.polygonActor = None

    def getNextNodeID(self):
        self.nextNodeID += 1
        return self.nextNodeID
    
    def createGrid(self, timeSlices, spatialPoints):
        if timeSlices < 3:
            print("ERROR: at least 3 timeslices needed!")

        if spatialPoints < 3:
            print("ERROR: at least 3 spatial point are needed!")

        self.numTimeSlices = timeSlices
        self.boundary = spatialPoints
        
        # Create the static grid
        timeList = []
        for ts in range(timeSlices):
            spaceList = []
            for ss in range(spatialPoints):
                n = self.createNewSpatialNode(ts)
                spaceList.append(n)
            timeList.append(spaceList)
            
        
        # Create spatial edges
        for ts in range(timeSlices):
            for ss in range(spatialPoints-1):
                timeList[ts][ss].right_neighbour = timeList[ts][ss+1]
            # timeList[ts][spatialPoints-1].right_neighbour = timeList[ts][0]  # For periodicity

        # Create time-like edges
        for ts in range(timeSlices-1):
            for ss in range(spatialPoints-1):
                cnode = timeList[ts][ss]
                f1node = timeList[ts+1][ss]
                f2node = timeList[ts+1][ss+1]
                
                cnode.future_nodes.append(f1node)
                f1node.past_nodes.append(cnode)
                
                cnode.future_nodes.append(f2node)
                f2node.past_nodes.append(cnode)
            
            # We have to handle the last node differently for periodicity
            cnode = timeList[ts][spatialPoints-1]
            f1node = timeList[ts+1][spatialPoints-1]
            # f2node = timeList[ts+1][0]      # For periodicity
        
            cnode.future_nodes.append(f1node)
            f1node.past_nodes.append(cnode)
                
            # cnode.future_nodes.append(f2node)      # For periodicity
            # f2node.past_nodes.insert(0, cnode)     # For periodicity
         
    def getVolumePerTimeSlice(self):
        '''
        Get the number of triangles = 2 * number of spatial nodes
        '''
        volumes = [0]*self.numTimeSlices
        
        for key, value in self.spatialNodes.items():
            volumes[value.t]+=1
            
        return volumes * 2
    
    def getTotalVolume(self):
        return (len(self.spatialNodes) - self.boundary) * 2
    
    def getRndNode(self):
        return random.choice(list(self.spatialNodes.values()))
    
    def stepUniverse(self):
        if random.random() > 0.5:
            if self.move():
                self.moveCounter += 1
        else:
            if self.inverseMove():
                self.invMoveCounter += 1
    
    def createNewSpatialNode(self, time):
        n = SpatialNode(self.getNextNodeID(), time)
        self.spatialNodes[n.ID] = n
        return n
    
    def eraseSpatialNode(self, node):
        self.spatialNodes.pop(node.ID)
    
    def move(self):
        # Get random node (not on boundary!)
        n = self.getRndNode()
        while n.t == 0 or n.t == (self.numTimeSlices-1) or n.right_neighbour == None:
            n = self.getRndNode()
        
        nf = len(n.future_nodes)
        np = len(n.past_nodes)
        
        # 1-based index for math
        splitf = random.randint(1, nf)
        splitp = random.randint(1, np)

        nf_new = nf - splitf + 1
        np_new = np - splitp + 1
        
        Nv = float(len(self.spatialNodes))
        
        Pmove = Nv / (Nv+1.0) * nf * np / (nf_new + np_new) * self.boltzmannWeight_neg
        
        if random.random() > Pmove:
            return False     
       
        # Create new node
        new_n = self.createNewSpatialNode(n.t)
        
        # Set spatial links
        new_n.right_neighbour = n.right_neighbour
        new_n.left_neighbour = n
        n.right_neighbour = new_n
        
        # Set timelike links
        new_n.future_nodes = n.future_nodes[splitf-1:]
        n.future_nodes = n.future_nodes[:splitf]
        
        new_n.past_nodes = n.past_nodes[splitp-1:]
        n.past_nodes = n.past_nodes[:splitp]
        
        # Set links for timelike neighbours
        #Insert split edges
        split_f_n = new_n.future_nodes[0]
        loc = split_f_n.past_nodes.index(n)
        split_f_n.past_nodes.insert(loc+1, new_n)
        
        split_p_n = new_n.past_nodes[0]
        loc = split_p_n.future_nodes.index(n)
        split_p_n.future_nodes.insert(loc+1, new_n)
        
        # Redirect moved old edges
        for f_ns in new_n.future_nodes[1:]:
            loc = f_ns.past_nodes.index(n)
            f_ns.past_nodes[loc] = new_n
            
        for p_ns in new_n.past_nodes[1:]:
            loc = p_ns.future_nodes.index(n)
            p_ns.future_nodes[loc] = new_n     
            
        return True
        
        
        
    def inverseMove(self):
        # Get random node (not on boundary!)
        n = self.getRndNode()
        while n.t == 0 or n.t == (self.numTimeSlices-1) or n.right_neighbour == None:
            n = self.getRndNode()

        old_n = n.right_neighbour
        
        nf_new = len(old_n.future_nodes)
        np_new = len(old_n.past_nodes)
        
        PinvMove = 1.0 / (nf_new + np_new) * self.boltzmannWeight_pos
       
        # Reject this sample
        if random.random() > PinvMove:
            return False          

        # Set spatial links
        n.right_neighbour = old_n.right_neighbour
        
        # Set timelike links
        n.future_nodes.extend(old_n.future_nodes[1:])
        n.past_nodes.extend(old_n.past_nodes[1:])
        
        
        #  Set links for timelike neighbours
        # Remove split edges
        split_f_n = old_n.future_nodes[0]
        split_f_n.past_nodes.remove(old_n)

        split_p_n = old_n.past_nodes[0]
        split_p_n.future_nodes.remove(old_n)
        
        # Redirect moved old edges
        for f_ns in old_n.future_nodes[1:]:
            loc = f_ns.past_nodes.index(old_n)
            f_ns.past_nodes[loc] = n
            
        for p_ns in old_n.past_nodes[1:]:
            loc = p_ns.future_nodes.index(old_n)
            p_ns.future_nodes[loc] = n     

        # Erase node from pool
        self.eraseSpatialNode(old_n)
        
        return True

    def listIDs(self, node_list):
        return [i.ID for i in node_list]

    def printUniverse(self):
        print("*********************")
        print("Times:", self.numTimeSlices)
        print("Original spatial points:", self.boundary)
        print("Volume:", self.getTotalVolume() )

    def saveUniverseCSV(self, fileName):
        # save nodes
        with open(fileName+"_nodes.csv", 'w') as fNodes:
            fNodes.write("Id,Label\n")
            for ID, node in self.spatialNodes.iteritems():
                fNodes.write(str(node.ID)+"," +str(node.t)+'\n')
                
        
        # save edges
        with open(fileName+"_edges.csv", 'w') as fNodes:
            fNodes.write("Source,Target,Type,Label\n")
            for ID, node in self.spatialNodes.iteritems():
                for e in node.past_nodes:
                    fNodes.write(str(node.ID)+"," +str(e.ID)+',Undirected,Time\n')
                if node.right_neighbour is not None:
                    fNodes.write(str(node.ID)+"," +str(node.right_neighbour.ID)+',Undirected,Space\n')
                

    def createVTKDataObject(self):
        import vtk
        pts = vtk.vtkPoints()
        lines = vtk.vtkCellArray()
        triangles = vtk.vtkCellArray()
        colors = vtk.vtkUnsignedCharArray()
        colors.SetNumberOfComponents(3)
        colors.SetName("Colors")

        #Get radius list
        numPerTime = self.getVolumePerTimeSlice()

        #Get a node on t=0
        pivot_node = None
        for ID, node in self.spatialNodes.items():
            if node.t == 0:
                pivot_node = node
                break
        
        # Create nodes
        for i in range(self.numTimeSlices):
            numPts = numPerTime[i]
            r = numPts / (2.0 * math.pi)
            x0 = 0.0
            y0 = 0.0
            
            pivot_space = pivot_node
            for j in range(numPts):
                x = x0 + r * math.cos(2 * math.pi * j / numPts)
                y = y0 + r * math.sin(2 * math.pi * j / numPts)
                pts.InsertPoint(pivot_space.ID, (x,y,float(pivot_space.t)))
                pivot_space = pivot_space.right_neighbour
            
            if i < self.numTimeSlices-1:
                pivot_node = pivot_node.future_nodes[0]
        
        
        # Add edges
        for ID, node in self.spatialNodes.items():
            if node.right_neighbour is not None:
                line = vtk.vtkLine()
                line.GetPointIds().SetId(0,node.ID) 
                line.GetPointIds().SetId(1,node.right_neighbour.ID)
                lines.InsertNextCell(line)
                colors.InsertNextTuple3(0,0,0)
            for e in node.future_nodes:
                line = vtk.vtkLine()
                line.GetPointIds().SetId(0,e.ID) 
                line.GetPointIds().SetId(1,node.ID)
                lines.InsertNextCell(line)
                colors.InsertNextTuple3(0,0,0)
                
        # Add triangles
        for ID, node in self.spatialNodes.items():
            for e in node.past_nodes:
                if e.right_neighbour is not None:
                    if node in e.right_neighbour.future_nodes:
                        triangle = vtk.vtkTriangle()
                        triangle.GetPointIds().SetId(0,node.ID)
                        triangle.GetPointIds().SetId(1,e.ID)
                        triangle.GetPointIds().SetId(2,e.right_neighbour.ID)
                        triangles.InsertNextCell(triangle)
                        colors.InsertNextTuple3(255,0,0)

            for e in node.future_nodes:
                if e.right_neighbour is not None:
                    if node in e.right_neighbour.past_nodes:
                        triangle = vtk.vtkTriangle()
                        triangle.GetPointIds().SetId(0,node.ID)
                        triangle.GetPointIds().SetId(1,e.ID)
                        triangle.GetPointIds().SetId(2,e.right_neighbour.ID)
                        triangles.InsertNextCell(triangle)
                        colors.InsertNextTuple3(0,0,255)
                
        
        # Create a polydata to store everything in
        linesPolyData = vtk.vtkPolyData()
 
        linesPolyData.SetPoints(pts)
        linesPolyData.SetLines(lines)
        linesPolyData.SetPolys(triangles)
        linesPolyData.GetCellData().SetScalars(colors)
        linesPolyData.Modified()
        
        return linesPolyData

    def createNewActor(self):
        polyData = self.createVTKDataObject()
        
        polygonMapper = vtk.vtkPolyDataMapper()

        polygonMapper.SetInputData(polyData)
        polygonMapper.Update()

        polygonActor = vtk.vtkActor()
        polygonActor.SetMapper(polygonMapper)

        return polygonActor

    def updateActor(self):
        self.render.RemoveActor(self.polygonActor)
        self.polygonActor = self.createNewActor()
        self.render.AddActor(self.polygonActor)

    def progressUniverse(self, steps):
        start = time.time()

        for i in range(steps): 
            self.stepUniverse()

        end = time.time()
        
        print("Progressing the Universe %d steps took %d seconds"%(steps, end-start))
        print("Total moves/Inv. moves:", u.moveCounter, "/", u.invMoveCounter)
        print("Total volume:", u.getTotalVolume() )

    def drawUniverseVTK(self):
        
        axesActor = vtk.vtkAxesActor()
        axes = vtk.vtkOrientationMarkerWidget()
        axes.SetOrientationMarker(axesActor)

        self.render = vtk.vtkRenderer()    

        # Get the starting version of the universe
        self.polygonActor = self.createNewActor()

        self.render.AddActor(self.polygonActor)
        self.render.SetBackground(0.1, 0.2, 0.4)
 
        self.render.ResetCamera()
 
        renWin = vtk.vtkRenderWindow()
        renWin.AddRenderer(self.render)
        renWin.SetSize(450, 350)
        renWin.SetWindowName("Causal Dynamical Triangulation 1+1D")
 
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
        iren.AddObserver("KeyPressEvent", self.keyPress)
        iren.SetRenderWindow(renWin)
        axes.SetInteractor(iren)
        axes.EnabledOn()
        axes.InteractiveOn()

        iren.Initialize()
        iren.Start()
    
    def keyPress(self,obj,event):

        key = obj.GetKeySym()
        #print(key)
        if key == "s":
            self.progressUniverse(progressSteps)
            self.updateActor()


    def saveUniverseVTK(self, fileName):
        
        polyData = self.createVTKDataObject()
        
        writer = vtk.vtkXMLPolyDataWriter();
        writer.SetFileName(fileName);
        writer.SetInputData(polyData)
        writer.Write()
        




if __name__ == "__main__":
    #random.seed(0)

    progressSteps = 50

    u = Universe()
    u.createGrid(32, 64)

    u.printUniverse()
    u.drawUniverseVTK()

    #progressUniverse(u, 500)

    #print("Runtime:", end-start, "s. It:", (end-start)/float(u.moveCounter+u.invMoveCounter) )

    #print (u.getVolumePerTimeSlice())
    #u.saveUniverseCSV("test1_universe") 
    #u.saveUniverseVTK("test1_universe.vtp")