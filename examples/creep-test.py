# to be able to add sofa objects you need to first load the plugins that implement them.
import SofaRuntime

# to add elements like Node or objects
import Sofa.Core
import math
import numpy as np
from scipy import signal

# import os
# path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'plot')
# if not os.path.isdir(path):
# 	if os.path.isfile(path):
# 		raise ValueError(f"path {path} already exist and is a file")
# 	else:
# 		os.mkdir(path)

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
		#self.tau = self.node.cylinder.FEM.ParameterSet.value[2] 
		epsilon = (self.node.cylinder.tetras.position.value[4][2]-self.lin)/self.lin

## IN THIS CODE WE WILL DO A CREEP TEST, SO WE WILL APPLY A STEP AS INPUT, USING THE CONSTANTFORCEFIELD LOAD.
 
##	STEP SIGNAL 
		self.node.cylinder.CFF.totalForce.value = [0,0,((3.141592653589793e2/2))*np.heaviside(self.time, 0.0)] ## amplitude in force calculated with Matlab --> Step function 
## Ramp

		#self.node.cylinder.CFF.totalForce.value =  [0,0,((3.141592653589793e2/2))*(1-abs(signal.sawtooth(2 * np.pi*7*self.time)))]## amplitude in force calculated with Matlab 
		#if(self.time>=0.1):
		#	self.node.cylinder.CFF.totalForce.value =  [0,0,((3.141592653589793e2/2))]## amplitude in force calculated with Matlab 



		print(epsilon*100)
		#file1 = open(os.path.join(path,"SLS_Maxwell_cyclic1.txt"),"a")
		#file1.write(str(self.time)+' '+str(self.node.cylinder.FEM.stressVonMisesElement.value[4])+' '+str(epsilon*100)+ '\n' )
		#file1.close()



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
	rootNode.addObject("RequiredPlugin", name="SofaViscoElastic")


	rootNode.gravity=[0,9.810,0]
	rootNode.dt = (1e6/(20e6*100))
	rootNode.name = 'rootNode'
	rootNode.addObject('DefaultAnimationLoop', computeBoundingBox="0")
	rootNode.addObject('BlockGaussSeidelConstraintSolver', tolerance=1e-24, maxIterations=1000)
	rootNode.addObject('OglSceneFrame', style='Arrows', alignment='TopRight')





## SLS-MAXWELL FIRST ORDER Cylinder: Material F5000- R0.5  (BLUE)
	cylinder = rootNode.addChild('cylinder')

	cylinder.addObject('EulerImplicitSolver', name="Solver",rayleighMass = 0.0, rayleighStiffness = 0.0, firstOrder = True, trapezoidalScheme = False)

	cylinder.addObject('CGLinearSolver', name="ItSolver", iterations="2500", tolerance="1e-30", threshold = '1e-12')
	cylinder.addObject('MeshVTKLoader', name='loader', filename='mesh/cylinder5296.vtk', translation = [0, 0.0, 0])
	cylinder.addObject('MechanicalObject', name='tetras', template='Vec3d', src = '@loader')
	cylinder.addObject('TetrahedronSetTopologyContainer', name="topo", src ='@loader')
	cylinder.addObject('TetrahedronSetTopologyModifier' ,  name="Modifier")
	cylinder.addObject('TetrahedronSetGeometryAlgorithms', template="Vec3d" ,name="GeomAlgo")

	cylinder.addObject('UniformMass', totalMass="1e-6", src = '@topo')

	cylinder.addObject('BoxROI', name='boxROI1',box="-0.011 -0.011 -0.001  0.011 0.011 0.001", drawBoxes=True)

	cylinder.addObject('FixedProjectiveConstraint', name = 'fix1', indices = '@boxROI1.indices')

	cylinder.addObject('BoxROI', name="boxToPull", box=[-0.011, -0.011, 0.05, 0.011, 0.011, 0.051], drawBoxes=False)

	G1 = 70e6
	G2 = 20e6
	G3 = 10e6
	tau1 = 1e6/G1
	tau2 = 1e6/G2
	tau3 = 1e6/G3
	nu = 0.44 ## Poisson's ratio
	lamb = 2*G1*nu/(1-2*nu) ## relationship from wikipedia: https://en.wikipedia.org/wiki/Lam%C3%A9_parameters

	cylinder.addObject('TetrahedronViscoelasticityFEMForceField', template='Vec3d', name='FEM', src ='@topo',materialName="SLSMaxwellFirstOrder", ParameterSet= str(G1)+' '+str(G2)+' '+str(tau2)+' '+str(lamb))
	
	cylinder.addObject('ConstantForceField', name = "CFF", listening = True, totalForce =[0,0,0],template="Vec3d", src= "@topo", indices = cylinder.boxToPull.indices.linkpath) 
	
	cylinder.addObject(CylinderController(node=rootNode, pos3 = rootNode.cylinder.tetras.position.value))


	modelVisu3 = cylinder.addChild('visu')
	modelVisu3.addObject('MeshSTLLoader', name='loader', filename='mesh/cylinder5296.stl', translation = [0.0,0.0,0.0])
	modelVisu3.addObject('OglModel', src='@loader', color=[1,1,1,1])
	modelVisu3.addObject('BarycentricMapping')







	return rootNode
