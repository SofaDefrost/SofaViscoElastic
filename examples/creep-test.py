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



		self.lin = self.node.cilinder.tetras.position.value[self.posmax1][2]
	
		file1 = open(path + "SLS_Maxwell.txt","w")
		file1.write(str(0.0)+' '+str(0.0)+ '\n' )
		file1.close()



	def onAnimateBeginEvent(self,event):
		#print(self.node.cilinder.FEM.getF(369))
		self.time = self.node.time.value
		epsilon = (self.node.cilinder.tetras.position.value[self.posmax1][2]-self.lin)/self.lin
		#self.node.cilinder.CFF.totalForce.value = [0,0,(3.141592653589793e4)*(signal.square(2 * np.pi * 10 * self.time))] ## --> Square wave function
		#self.node.cilinder.CFF.totalForce.value = [0,0,(3.141592653589793e4)*np.sin(2*np.pi*0.6*self.time)] ## --> Sinusoidal function
		self.node.cilinder.CFF.totalForce.value = [0,0,3.141592653589793e4] ## Calculated with Matlab --> Step function 
		print(epsilon*100)
		file1 = open(path + "SLS_Maxwell.txt","a")
		file1.write(str(self.time)+' '+str(self.node.cilinder.CFF.totalForce.value[2])+' '+str(epsilon*100)+ '\n' )
		file1.close()



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
	rootNode.addObject("RequiredPlugin", name="Sofa.Component.SolidMechanics.FEM.ViscoElastic")


	rootNode.gravity=[0,0,-9.81]
	rootNode.dt = (90/1400)
	rootNode.name = 'rootNode'
	rootNode.addObject('DefaultAnimationLoop', computeBoundingBox="0")
	rootNode.addObject('GenericConstraintSolver', tolerance=1e-12, maxIterations=10000)
	rootNode.addObject('OglSceneFrame', style='Arrows', alignment='TopRight')





## SLS-MAXWELL FIRST ORDER Cylinder: Material F5000- R0.5  (BLUE)
	cilinder = rootNode.addChild('cilinder')

	cilinder.addObject('EulerImplicitSolver', name="Solver",rayleighMass = 0.0, rayleighStiffness = 0.0, firstOrder = False)

	cilinder.addObject('CGLinearSolver', name="ItSolver", iterations="2500", tolerance="1e-15", threshold = '1e-30')
	cilinder.addObject('MeshVTKLoader', name='loader', filename='mesh/cylinder1513.vtk', translation = [0, 0.0, 0])
	cilinder.addObject('MechanicalObject', name='tetras', template='Vec3d', src = '@loader')
	cilinder.addObject('TetrahedronSetTopologyContainer', name="topo", src ='@loader')
	cilinder.addObject('TetrahedronSetTopologyModifier' ,  name="Modifier")
	cilinder.addObject('TetrahedronSetGeometryAlgorithms', template="Vec3d" ,name="GeomAlgo")

	cilinder.addObject('UniformMass', totalMass="0.013", src = '@topo')

	cilinder.addObject('BoxROI', name='boxROI1',box="-0.011 -0.011 -0.001  0.011 0.011 0.001", drawBoxes=True)
	cilinder.addObject('FixedConstraint', indices = '@boxROI1.indices')
	E0 = 70e9
	E1 = 20e9
	tau1 = (90)/(70*20)
	nu = 0.44
	cilinder.addObject('TetrahedronViscoelasticityFEMForceField', template='Vec3d', name='FEM', src ='@topo',materialName="SLSMaxwellFirstOrder", ParameterSet=str(E0)+' '+str(E1)+' '+str(tau1)+' '+str(nu))

	cilinder.addObject('ConstantForceField', name = "CFF", listening = True, totalForce =[0,0,0],template="Vec3d", src= "@topo", indices =" 4 5 6 7 20 21 22 23 24 25 26 27 28 29 30 31 48 49 50 67 68 69 86 87 88 105 106 107 143 144 146 149 160 161 178 186 190 202 203 205 208 219 220 237 245 249 261 262 264 267 278 279 296 304 308 320 321 323 326 337 338 355 363 367 369 370 371 372 373 374 375 376 377 378 379 380 381 382 383 384 385 386 387 388 389 390 391 392 393 396 404 408 409 411 412 414 415 421 455 458 459 460 461 468 469 470 471 476 477 478 479 483 ") 
	cilinder.addObject(CylinderController(node=rootNode, pos3 = rootNode.cilinder.tetras.position.value))


	modelVisu3 = cilinder.addChild('visu')
	modelVisu3.addObject('MeshSTLLoader', name='loader', filename='mesh/cylinder5296.stl', translation = [0.0,0.0,0.0])
	modelVisu3.addObject('OglModel', src='@loader', color=[0,0,1,1])
	modelVisu3.addObject('BarycentricMapping')







	return rootNode
