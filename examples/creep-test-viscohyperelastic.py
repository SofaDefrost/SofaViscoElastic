# to be able to add sofa objects you need to first load the plugins that implement them.
import SofaRuntime

# to add elements like Node or objects
import Sofa.Core
import math
import numpy as np
from scipy import signal



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
		#self.node.cylinder.CFF.totalForce.value = [0, 0, ((3.141592653589793e3))*7.8*np.heaviside(self.time, 0.0)] ## amplitude in force calculated with Matlab --> Step function 

		#if(self.time >= 1):
		#	self.node.cylinder.CFF.totalForce.value = [0, 0, 0] ## amplitude in force calculated with Matlab --> Step function 

		#self.node.cylinder.CFF.totalForce.value = [0, 0, 0] ## amplitude in force calculated with Matlab --> Step function 

		
##	AMPLITUDE SWEEP 
		self.node.cylinder.CFF.totalForce.value = [0,0,((3.141592653589793e2/2))*1*np.heaviside(self.time, 0.0)] ## amplitude in force calculated with Matlab --> Step function 

		##if(self.time >= self.tau*100 and self.time<= self.tau*150):
		##	self.node.cylinder.CFF.totalForce.value = [0, 0, 0] ## amplitude in force calculated with Matlab --> Step function
		

		##if(self.time >= self.tau*150 and self.time <= self.tau*200):		
		##	self.node.cylinder.CFF.totalForce.value = [0,0,((3.141592653589793e3/2))*3*np.heaviside(self.time, 0.0)] ## amplitude in force calculated with Matlab --> Step function 
		
		##if(self.time >= self.tau*200 and self.time <= self.tau*250):		
		#	self.node.cylinder.CFF.totalForce.value = [0,0,0] ## amplitude in force calculated with Matlab --> Step function 

		#if(self.time >= self.tau*250 and self.time <= self.tau*300):		
		#	self.node.cylinder.CFF.totalForce.value = [0,0,((3.141592653589793e3/2))*5*np.heaviside(self.time, 0.0)] ## amplitude in force calculated with Matlab --> Step function 
		
		#if(self.time >= self.tau*300 and self.time <= self.tau*350):		
		#	self.node.cylinder.CFF.totalForce.value = [0,0,0] ## amplitude in force calculated with Matlab --> Step function 		

		#if(self.time >= self.tau*350 and self.time <= self.tau*400):		
		#	self.node.cylinder.CFF.totalForce.value = [0,0,((3.141592653589793e3/2))*7*np.heaviside(self.time, 0.0)] ## amplitude in force calculated with Matlab --> Step function 
		
		#if(self.time >= self.tau*400):		
		#	self.node.cylinder.CFF.totalForce.value = [0,0,0] ## amplitude in force calculated with Matlab --> Step function 

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
	rootNode.addObject("RequiredPlugin", name="SofaViscoElastic")


	rootNode.gravity=[0,0,-9.81]
	rootNode.dt = (1e6/(10e6*100))
	rootNode.name = 'rootNode'
	rootNode.addObject('DefaultAnimationLoop', computeBoundingBox="0")
	rootNode.addObject('BlockGaussSeidelConstraintSolver', tolerance=1e-24, maxIterations=100000000)
	rootNode.addObject('OglSceneFrame', style='Arrows', alignment='TopRight')





## SLS-MAXWELL FIRST ORDER Cylinder: Material F5000- R0.5  (BLUE)
	cylinder = rootNode.addChild('cylinder')

	cylinder.addObject('EulerImplicitSolver', name="Solver",rayleighMass = 0.0, rayleighStiffness = 0.0)

	cylinder.addObject('CGLinearSolver', name="ItSolver", iterations="250", tolerance="1e-12", threshold = '1e-12')
	cylinder.addObject('MeshVTKLoader', name='loader', filename='mesh/cylinder5296.vtk', translation = [0, 0.0, 0])
	cylinder.addObject('MechanicalObject', name='tetras', template='Vec3d', src = '@loader')
	cylinder.addObject('TetrahedronSetTopologyContainer', name="topo", src ='@loader')
	cylinder.addObject('TetrahedronSetTopologyModifier' ,  name="Modifier")
	cylinder.addObject('TetrahedronSetGeometryAlgorithms', template="Vec3d" ,name="GeomAlgo")

	cylinder.addObject('UniformMass', totalMass="1", src = '@topo')

	cylinder.addObject('BoxROI', name='boxROI1',box="-0.011 -0.011 -0.001  0.011 0.011 0.001", drawBoxes=True)
	cylinder.addObject('FixedProjectiveConstraint', indices = '@boxROI1.indices')
	cylinder.addObject('BoxROI', name="boxToPull", box=[-0.011, -0.011, 0.1, 0.011, 0.011, 0.101], drawBoxes=False)

	mu = 70e6
	E2 = 20e6
	E3 = 10e6
	alpha = 1
	tau1 = 1e6/mu
	tau2 = 1e6/E2
	tau3 = 1e6/E3
	k = 44e6
	cylinder.addObject('TetrahedronViscoHyperelasticityFEMForceField', template='Vec3d', name='FEM', src ='@topo',materialName="SLSNeoHookeanFirstOrder", ParameterSet= str(mu)+' '+str(mu)+' '+str(tau1)+' '+str(k))
	
	cylinder.addObject('ConstantForceField', name = "CFF", listening = True, totalForce =[0,0,0],template="Vec3d", src= "@topo", indices = cylinder.boxToPull.indices.linkpath) 
	
	cylinder.addObject(CylinderController(node=rootNode, pos3 = rootNode.cylinder.tetras.position.value))


	modelVisu3 = cylinder.addChild('visu')
	modelVisu3.addObject('MeshSTLLoader', name='loader', filename='mesh/cylinder5296.stl', translation = [0.0,0.0,0.0])
	modelVisu3.addObject('OglModel', src='@loader', color=[1,0,0,1])
	modelVisu3.addObject('BarycentricMapping')







	return rootNode
