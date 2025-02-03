import SofaRuntime
SofaRuntime.importPlugin("SofaComponentAll")

# to add elements like Node or objects
import Sofa.Core
root = Sofa.Core.Node()
import math 
import numpy as np

import os
path = os.path.dirname(os.path.abspath(__file__))+'/plot/'


class CylinderController(Sofa.Core.Controller):

	def __init__(self, *args, **kwargs):
		Sofa.Core.Controller.__init__(self,*args, **kwargs) #needed
		self.time = 0.0
		self.node = kwargs['node']
		self.pos3 = kwargs['pos']

		self.max1 = 0
		self.posmax1 = 0
		self.max2 = 0
		self.posmax2 = 0
		for j in range(0, len(self.pos3)):
			if self.pos3[j][2] >= self.max1 :
				self.max1 = self.pos3[j][2]
				self.posmax1 = j

		print(self.pos3[self.posmax1][2], self.posmax1)

		self.lin = self.node.cylinder.tetras.position.value[self.posmax1][2]



	def onAnimateBeginEvent(self,event):

	
		stress = (self.node.cylinder.FEM.CauchyStress.value)/1e6 ## MPa
		stressPerNode = np.array([[0,0,0,0],[0,0,0,0]])


				
		self.time = self.node.time.value
		epsilon = (self.node.cylinder.tetras.position.value[self.posmax1][2]-self.lin)/self.lin
		self.node.cylinder.visu.display.pointData.value = (stress[0:len(self.pos3),2])
		self.node.cylinder.visu.Map.min.value = self.node.cylinder.visu.display.currentMin.value
		self.node.cylinder.visu.Map.max.value = self.node.cylinder.visu.display.currentMax.value




## IN THIS CODE WE WILL DO A STRESS RELAXATION  TEST, SO WE WILL APPLY A STEP AS INPUT, USING THE POSITIONCONSTRAINT.



def createScene(rootNode):
	rootNode.addObject('RequiredPlugin', name='SofaPython3')
	rootNode.addObject('RequiredPlugin', name="Sofa.Component.Engine.Select")
	rootNode.addObject('RequiredPlugin', name="Sofa.Component.IO.Mesh")
	rootNode.addObject("RequiredPlugin", name="Sofa.Component.LinearSolver.Iterative")
	rootNode.addObject("RequiredPlugin", name="Sofa.Component.Mapping.Linear")
	rootNode.addObject("RequiredPlugin", name="Sofa.Component.Mass")
	rootNode.addObject("RequiredPlugin", name="Sofa.Component.ODESolver.Forward")
	rootNode.addObject("RequiredPlugin", name="Sofa.Component.Setting")
	rootNode.addObject("RequiredPlugin", name="Sofa.Component.SolidMechanics.FEM.HyperElastic")
	rootNode.addObject("RequiredPlugin", name="Sofa.Component.SolidMechanics.FEM.Elastic")	
	rootNode.addObject("RequiredPlugin", name="Sofa.Component.SolidMechanics.Spring")
	rootNode.addObject("RequiredPlugin", name="Sofa.Component.StateContainer")
	rootNode.addObject("RequiredPlugin", name="Sofa.Component.Topology.Container.Dynamic")
	rootNode.addObject("RequiredPlugin", name="Sofa.Component.Visual")
	rootNode.addObject("RequiredPlugin", name= "Sofa.GL.Component.Rendering3D")
	rootNode.addObject("RequiredPlugin", name="Sofa.Component.AnimationLoop")
	rootNode.addObject("RequiredPlugin", name="Sofa.Component.Constraint.Lagrangian.Solver")
	rootNode.addObject("RequiredPlugin", name="Sofa.Component.MechanicalLoad")
	rootNode.addObject("RequiredPlugin", name="Sofa.Component.LinearSolver.Direct")
	rootNode.addObject("RequiredPlugin", name="Sofa.Component.Constraint.Lagrangian.Correction")
	rootNode.addObject("RequiredPlugin", name = "Sofa.Component.Constraint.Projective")
	rootNode.addObject("RequiredPlugin", name="Sofa.Component.ODESolver.Backward")
	rootNode.addObject('RequiredPlugin', name='Sofa.GL.Component.Rendering2D') # Needed to use components [OglColorMap]
	rootNode.addObject('RequiredPlugin', name='SoftRobots') # Needed to use components [PositionConstraint] 

	rootNode.addObject('FreeMotionAnimationLoop')
	rootNode.addObject('GenericConstraintSolver', maxIterations=1e4, tolerance=1e-50)
	rootNode.gravity = [0,0,-9.81]
	rootNode.dt = (1e9/(20e9*1000))

	rootNode.addObject('VisualStyle', displayFlags='hideForceFields')
	rootNode.addObject('OglSceneFrame', style='Arrows', alignment='TopRight')	

	cylinder = rootNode.addChild('cylinder')

	cylinder.addObject('EulerImplicitSolver', name="Solver",rayleighMass = 0.0, rayleighStiffness = 0.0, firstOrder = True, trapezoidalScheme = False)
	cylinder.addObject('SparseLDLSolver', name="directsolver")

	cylinder.addObject('MeshVTKLoader', name='loader', filename='mesh/cylinder9178.vtk',)
	cylinder.addObject('MechanicalObject', name='tetras', template='Vec3d', src = '@loader', rotation = [0, 0 , 0])
	cylinder.addObject('TetrahedronSetTopologyContainer', name="topo", src ='@loader')
	cylinder.addObject('TetrahedronSetTopologyModifier' ,  name="Modifier")
	cylinder.addObject('TetrahedronSetGeometryAlgorithms', template="Vec3d" ,name="GeomAlgo")

	cylinder.addObject('UniformMass', totalMass="0.0126", src = '@topo')
	E1 = 70e9
	tau1 = 1e9/E1
	E2 = 20e9
	tau2 = 1e9/E2
	E3 = 10e9
	tau3 = 1e9/E3
	nu = 0.44

	cylinder.addObject('TetrahedronViscoelasticityFEMForceField', template='Vec3d', name='FEM', src ='@topo',materialName="SLSKelvinVoigtSecondOrder", ParameterSet= str(E1)+' '+str(E2)+' '+str(tau2)+' '+str(E3)+' '+str(tau3)+' '+str(nu))

	cylinder.addObject('BoxROI', name='boxROI',box="-0.011 -0.011 -0.001  0.011 0.011 0.001", drawBoxes=True)
	cylinder.addObject('FixedProjectiveConstraint', indices = '@boxROI.indices')
	cylinder.addObject('BoxROI', name="boxToPull", box=[-0.011, -0.011, 0.1, 0.011, 0.011, 0.101], drawBoxes=True)
	cylinder.addObject('PartialFixedProjectiveConstraint', indices=cylinder.boxToPull.indices.linkpath, fixedDirections=[1, 1, 0])

##	STEP SIGNAL 	
	cylinder.addObject('PositionConstraint', name = 'displacement', indices=cylinder.boxToPull.indices.linkpath,valueType="displacement", value = 1e-4 , useDirections=[0, 0, 1])
	cylinder.addObject('LinearSolverConstraintCorrection')

	cylinder.addObject(CylinderController(node=rootNode, pos = rootNode.cylinder.tetras.position))



	modelVisu1 = cylinder.addChild('visu')
	modelVisu1.addObject('DataDisplay', name = 'display') 
	modelVisu1.addObject('OglColorMap', name = 'Map', colorScheme = "Blue to Red", showLegend = True )
	modelVisu1.addObject('IdentityMapping',input="@..", output="@.")
	return rootNode