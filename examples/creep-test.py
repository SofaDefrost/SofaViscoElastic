# to be able to add sofa objects you need to first load the plugins that implement them.
# For simplicity you can load the plugin "SofaComponentAll" that will load all most
# common sofa objects.
import SofaRuntime
SofaRuntime.importPlugin("SofaComponentAll")

# to add elements like Node or objects
import Sofa.Core
root = Sofa.Core.Node()
import SofaViscoElastic
import math 
import numpy as np
from scipy import signal

import os
path = os.path.dirname(os.path.abspath(__file__))+'/plot/'

class CylinderController(Sofa.Core.Controller):

	def __init__(self, *args, **kwargs):
		Sofa.Core.Controller.__init__(self,*args, **kwargs) #needed
		self.time = 0.0
		self.node = kwargs['node']
		self.pos3 = kwargs['pos3']

		self.max1 = 0
		self.posmax1 = 0
		self.max2 = 0
		self.posmax2 = 0

		for j in range(0, len(self.pos3)):
			if self.pos3[j][2] >= self.max1:
				self.max1 = self.pos3[j][2]
				self.posmax1 = j

		print(self.posmax1)

		self.lin = self.node.cylinder.tetras.position.value[self.posmax1][2]
	



	def onAnimateBeginEvent(self,event):
		self.time = self.node.time.value
		self.tau = self.node.cylinder.FEM.ParameterSet.value[2] 
		epsilon = (self.node.cylinder.tetras.position.value[4][2]-self.lin)/self.lin

## IN THIS CODE WE WILL DO A CREEP TEST, SO WE WILL APPLY A STEP AS INPUT, USING THE CONSTANTFORCEFIELD LOAD.
 

##	STEP SIGNAL 
		self.node.cylinder.CFF.totalForce.value = [0,0,3.141592653589793e4*np.heaviside(self.time, 0.0)] ## amplitude in force calculated with Matlab --> Step function 



		print(epsilon*100)



def createScene(rootNode):

	rootNode.addObject('VisualStyle', displayFlags='showVisualModels hideBehaviorModels hideCollisionModels hideBoundingCollisionModels hideForceFields hideInteractionForceFields hideWireframe')
	
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
	rootNode.addObject("RequiredPlugin", name="Sofa.GL.Component.Rendering3D")
	rootNode.addObject("RequiredPlugin", name="Sofa.Component.AnimationLoop")
	rootNode.addObject("RequiredPlugin", name="Sofa.Component.Constraint.Lagrangian.Solver")
	rootNode.addObject("RequiredPlugin", name="Sofa.Component.MechanicalLoad")
	rootNode.addObject("RequiredPlugin", name="Sofa.Component.LinearSolver.Direct")
	rootNode.addObject("RequiredPlugin", name="Sofa.Component.Constraint.Lagrangian.Correction")
	rootNode.addObject("RequiredPlugin", name = "Sofa.Component.Constraint.Projective")
	rootNode.addObject("RequiredPlugin", name="Sofa.Component.ODESolver.Backward")


	rootNode.gravity=[0,0,-9.81]
	rootNode.dt = (1e9/(20e9*100))
	rootNode.name = 'rootNode'
	rootNode.addObject('DefaultAnimationLoop', computeBoundingBox="0")
	rootNode.addObject('GenericConstraintSolver', tolerance=1e-24, maxIterations=100000000)
	rootNode.addObject('OglSceneFrame', style='Arrows', alignment='TopRight')





## SLS-MAXWELL FIRST ORDER Cylinder: Material F5000- R0.5  (BLUE)
	cylinder = rootNode.addChild('cylinder')

	cylinder.addObject('EulerImplicitSolver', name="Solver",rayleighMass = 0.0, rayleighStiffness = 0.0, firstOrder = True, trapezoidalScheme = False)

	cylinder.addObject('CGLinearSolver', name="ItSolver", iterations="25000000", tolerance="1e-15", threshold = '1e-30')
	cylinder.addObject('MeshVTKLoader', name='loader', filename='mesh/cylinder5296.vtk', translation = [0, 0.0, 0])
	cylinder.addObject('MechanicalObject', name='tetras', template='Vec3d', src = '@loader')
	cylinder.addObject('TetrahedronSetTopologyContainer', name="topo", src ='@loader')
	cylinder.addObject('TetrahedronSetTopologyModifier' ,  name="Modifier")
	cylinder.addObject('TetrahedronSetGeometryAlgorithms', template="Vec3d" ,name="GeomAlgo")

	cylinder.addObject('UniformMass', totalMass="1e-6", src = '@topo')

	cylinder.addObject('BoxROI', name='boxROI1',box="-0.011 -0.011 -0.001  0.011 0.011 0.001", drawBoxes=True)
	cylinder.addObject('FixedConstraint', indices = '@boxROI1.indices')
	cylinder.addObject('BoxROI', name="boxToPull", box=[-0.011, -0.011, 0.1, 0.011, 0.011, 0.101], drawBoxes=True)
	cylinder.addObject('PartialFixedConstraint', indices=cylinder.boxToPull.indices.linkpath, fixedDirections=[1, 1, 0])

	E1 = 70e9
	E2 = 20e9
	tau1 = 1e9/E1
	tau2 = 1e9/E2
	nu = 0.44
	cylinder.addObject('TetrahedronViscoelasticityFEMForceField', template='Vec3d', name='FEM', src ='@topo',materialName="Burgers", ParameterSet= str(E1)+' '+str(tau1)+' '+str(E2)+' '+str(tau2)+' '+str(nu))
	
	cylinder.addObject('ConstantForceField', name = "CFF", listening = True, totalForce =[0,0,0],template="Vec3d", src= "@topo", indices = cylinder.boxToPull.indices.linkpath) 
	
	cylinder.addObject(CylinderController(node=rootNode, pos3 = rootNode.cylinder.tetras.position.value))


	modelVisu3 = cylinder.addChild('visu')
	modelVisu3.addObject('MeshSTLLoader', name='loader', filename='mesh/cylinder5296.stl', translation = [0.0,0.0,0.0])
	modelVisu3.addObject('OglModel', src='@loader', color=[0,0,1,1])
	modelVisu3.addObject('BarycentricMapping')







	return rootNode
