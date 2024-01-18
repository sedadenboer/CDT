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
    """Represents a spatial node in the universe, corresponding to a vertex in the triangulation.

    Attributes:
        ID: Unique identifier for the node.
        t: Time slice to which the node belongs.
        future_nodes: List of nodes in the future (next time slice) connected to this node.
        past_nodes: List of nodes in the past (previous time slice) connected to this node.
        right_neighbour: The spatially adjacent node to the right (connected in space).
        boundary: Indicates whether the node is on the boundary (0 or 1).
    """

    def __init__(self, ID, time):
        self.ID = ID
        self.t = time
        self.future_nodes = []
        self.past_nodes = []
        self.right_neighbour = None
        self.left_neighbour = None 
        self.boundary = 0  

class Universe(object):
    """Represents the entire universe and manages the collection of spatial nodes.

    Attributes:
        spatialNodes (dict): Dictionary to store spatial nodes with their IDs as keys.
        numTimeSlices (int): Number of time slices in the universe.
        nextNodeID (int): Counter for generating unique IDs for nodes.
        scaledCosmologicalConstant (float): Constant used in Boltzmann weights for Metropolis moves.
        boltzmannWeight_neg (float): Boltzmann weight for negative moves.
        boltzmannWeight_pos (float): Boltzmann weight for positive moves.
        moveCounter (int): Counter for positive Metropolis moves.
        invMoveCounter (int): Counter for inverse Metropolis moves.
        render (vtkRenderer): VTK renderer for visualization.
        polygonActor (vtkActor): VTK actor for visualizing the triangulation.
    """

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
        """Increments and returns the next available node ID."""
        self.nextNodeID += 1
        return self.nextNodeID
    
    def createGrid(self, timeSlices, spatialPoints):
        """Creates the initial spatial grid for the universe.

        Args:
            timeSlices (int): Number of time slices in the universe.
            spatialPoints (int): Number of spatial points (vertices) in each time slice.

        Raises:
            ValueError: If the number of time slices or spatial points is less than 3.
        """
        if timeSlices < 3:
            print("ERROR: at least 3 timeslices needed!")

        if spatialPoints < 3:
            print("ERROR: at least 3 spatial point are needed!")

        self.numTimeSlices = timeSlices
        self.boundary = spatialPoints
        
        # Create the static grid
        timeList = []
        # Create spatial nodes
        for ts in range(timeSlices):
            spaceList = []
            # Create spatial nodes
            for ss in range(spatialPoints):
                # Create new node
                n = self.createNewSpatialNode(ts)
                spaceList.append(n)
            timeList.append(spaceList)
            
        
        # Create spatial edges
        for ts in range(timeSlices):
            for ss in range(spatialPoints-1):
                # Set spatial links
                timeList[ts][ss].right_neighbour = timeList[ts][ss+1]
            # timeList[ts][spatialPoints-1].right_neighbour = timeList[ts][0]  # For periodicity

        # Create time-like edges
        for ts in range(timeSlices-1):
            for ss in range(spatialPoints-1):
                # Set timelike links
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
        """Computes the number of triangles (2 * number of spatial nodes) per time slice.

        Returns:
            list: A list representing the number of triangles per time slice.
        """
        volumes = [0]*self.numTimeSlices
        
        for key, value in self.spatialNodes.items():
            volumes[value.t]+=1
            
        return volumes * 2
    
    def getTotalVolume(self):
        """Computes the total volume of the universe."""
        return (len(self.spatialNodes) - self.boundary) * 2
    
    def getRndNode(self):
        """Returns a random spatial node from the universe.
        
        Returns:
            SpatialNode: The randomly selected spatial node.
        """
        return random.choice(list(self.spatialNodes.values()))
    
    def stepUniverse(self):
        """Performs a single Metropolis step on the universe."""
        if random.random() > 0.5:
            if self.move():
                self.moveCounter += 1
        else:
            if self.inverseMove():
                self.invMoveCounter += 1
    
    def createNewSpatialNode(self, time):
        """Creates a new spatial node and adds it to the universe.

        Args:
            time (int): Time slice to which the node belongs.

        Returns:
            SpatialNode: The newly created spatial node.
        """
        n = SpatialNode(self.getNextNodeID(), time)
        self.spatialNodes[n.ID] = n
        return n
    
    def eraseSpatialNode(self, node):
        """Removes the given spatial node from the universe."""
        self.spatialNodes.pop(node.ID)
    
    def move(self):
        """Implements a Metropolis move by splitting and reconnecting nodes.

        Returns:
            bool: True if the move was accepted, False otherwise.
        """
        # Get random node (not on boundary!)
        n = self.getRndNode()
        while n.t == 0 or n.t == (self.numTimeSlices-1) or n.right_neighbour == None:
            n = self.getRndNode()
        
        # Get number of future and past nodes
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
        """Implements a Metropolis move by merging and reconnecting nodes.
        
        Returns:
            bool: True if the move was accepted, False otherwise.
        """
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
        """Returns a list of IDs for the given list of nodes.

        Args:
            node_list (list): List of nodes.

        Returns:
            list: List of IDs.
        """
        return [i.ID for i in node_list]

    def printUniverse(self):
        """Prints the universe to the console."""
        print("*********************")
        print("Times:", self.numTimeSlices)
        print("Original spatial points:", self.boundary)
        print("Volume:", self.getTotalVolume() )

    def saveUniverseCSV(self, fileName):
        """Saves the universe to a CSV file.

        Args:
            fileName (str): Name of the file to save to.
        
        Note:
            The saved CSV file can be imported into Gephi for visualization.
        """
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
        """Creates a VTK data object from the universe.

        Returns:
            vtkPolyData: The VTK data object representing the universe.
        """
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
        """Creates a VTK actor from the universe.

        Returns:
            vtkActor: The VTK actor representing the universe.
        """
        polyData = self.createVTKDataObject()
        
        polygonMapper = vtk.vtkPolyDataMapper()

        polygonMapper.SetInputData(polyData)
        polygonMapper.Update()

        polygonActor = vtk.vtkActor()
        polygonActor.SetMapper(polygonMapper)

        return polygonActor

    def updateActor(self):
        """Updates the VTK actor representing the universe."""
        self.render.RemoveActor(self.polygonActor)
        self.polygonActor = self.createNewActor()
        self.render.AddActor(self.polygonActor)

    def progressUniverse(self, steps):
        """Progresses the universe by the given number of steps.

        Args:
            steps (int): Number of steps to progress the universe.
        """
        start = time.time()

        for i in range(steps): 
            self.stepUniverse()

        end = time.time()
        
        print("Progressing the Universe %d steps took %d seconds"%(steps, end-start))
        print("Total moves/Inv. moves:", u.moveCounter, "/", u.invMoveCounter)
        print("Total volume:", u.getTotalVolume() )

    def drawUniverseVTK(self):
        """Draws the universe using VTK."""
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
        """Handles key presses in the VTK window."""	
        key = obj.GetKeySym()
        #print(key)
        if key == "s":
            self.progressUniverse(progressSteps)
            self.updateActor()


    def saveUniverseVTK(self, fileName):
        """Saves the universe to a VTK file."""
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