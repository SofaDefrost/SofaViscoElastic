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

import os
path = os.path.dirname(os.path.abspath(__file__))+'/plot/'

class CylinderController(Sofa.Core.Controller):

	def __init__(self, *args, **kwargs):
		Sofa.Core.Controller.__init__(self,*args, **kwargs) #needed
		self.time = 0.0
		self.node = kwargs['node']
		##self.Positions =  kwargs['position']

		self.CFF1 = kwargs['CFF1']
		self.CFF2 = kwargs['CFF2']
		self.CFF3 = kwargs['CFF3']
		self.CFF4 = kwargs['CFF4']

		self.pos1 = kwargs['pos1']
		self.pos2 = kwargs['pos2']
		self.pos3 = kwargs['pos3']
		self.pos4 = kwargs['pos4']

		self.max1 = 0
		self.posmax1 = 0
		self.max2 = 0
		self.posmax2 = 0

		for j in range(0, len(self.pos3)):
			if self.pos3[j][2] >= self.max1:
				self.max1 = self.pos3[j][2]
				self.posmax1 = j

		for k in range(0, len(self.pos4)):
			if self.pos4[j][2] >= self.max2:
				self.max2 = self.pos4[j][2]
				self.posmax2 = k


		#file1 = open("/home/pasquale/Script_Sofa/example/plot/posMaxwell.txt","w")
		#file1.write(str(0.0)+' '+str(self.pos3[self.posmax1][2])+ '\n' )
		#file1.close()

		#file2 = open("/home/pasquale/Script_Sofa/example/plot/posKelvin.txt","w")
		#file2.write(str(0.0)+' '+str(self.pos4[self.posmax2][2])+ '\n' )
		#file2.close()

	def onAnimateBeginEvent(self,event):
		print(self.node.cilinder3.FEM.getF(369))

		self.time += self.node.dt.value 
		f = 3; ## frequency Hz
		self.CFF1.force.value = [0,0,(0.1/2)-(0.1/2)*math.cos(2*math.pi*self.time*f)]
		self.CFF2.force.value = [0,0,(1/2)-(1/2)*math.cos(2*math.pi*self.time*f)]
		self.CFF3.force.value = [0,0,(1/2)-(1/2)*math.cos(2*math.pi*self.time*f)]
		self.CFF4.force.value = [0,0,(1/2)-(1/2)*math.cos(2*math.pi*self.time*f)]

		#file1 = open("/home/pasquale/Script_Sofa/example/plot/posMaxwell.txt","a")
		#file1.write(str(self.time)+' '+str(self.pos3[self.posmax1][2])+ '\n' )
		#file1.close()

		#file2 = open("/home/pasquale/Script_Sofa/example/plot/posKelvin.txt","a")
		#file2.write(str(self.time)+' '+str(self.pos4[self.posmax2][2])+ '\n' )
		#file2.close()



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


	rootNode.gravity=[0,0,0]
	rootNode.dt = 1e-2
	rootNode.name = 'rootNode'
	rootNode.addObject('DefaultAnimationLoop', computeBoundingBox="0")
	rootNode.addObject('GenericConstraintSolver', tolerance=1e-12, maxIterations=10000)


## MAXWELL FIRST ORDER Cylinder: Material F5000 R 0.5 (RED)

	rootNode.addObject('OglSceneFrame', style='Arrows', alignment='TopRight')

	cilinder1 = rootNode.addChild('cilinder1')

	cilinder1.addObject('EulerImplicitSolver', name="Solver")

	cilinder1.addObject('CGLinearSolver', name="ItSolver", iterations="250", tolerance="1e-15", threshold = '1e-40')

	cilinder1.addObject('MeshVTKLoader', name='loader', filename='mesh/cylinder1513.vtk',)
	cilinder1.addObject('MechanicalObject', name='tetras', template='Vec3d', src = '@loader', rotation = [0, 0 , 0])
	cilinder1.addObject('TetrahedronSetTopologyContainer', name="topo", src ='@loader')
	cilinder1.addObject('TetrahedronSetTopologyModifier' ,  name="Modifier")
	cilinder1.addObject('TetrahedronSetGeometryAlgorithms', template="Vec3d" ,name="GeomAlgo")

	cilinder1.addObject('UniformMass', totalMass="0.0126", src = '@topo')

	cilinder1.addObject('BoxROI', name='boxROI',box="-0.011 -0.011 -0.001  0.011 0.011 0.001", drawBoxes=True)
	cilinder1.addObject('FixedConstraint', indices = '@boxROI.indices')
	E1 = 1068000
	tau1 = 0.897
	k0 = 1/1.150e-3
	cilinder1.addObject('TetrahedronViscoelasticityFEMForceField', template='Vec3d', name='FEM', src ='@topo',materialName="MaxwellFirstOrder", ParameterSet= str(E1)+' '+str(tau1)+' '+str(k0))


	cilinder1.addObject('ConstantForceField', name = "CFF", listening = True, force =[0,0,0],indices ="4 5 6 7 20 21 22 23 28 29 30 31 68 69 87 88 106 107 143 144 146 149 203 205 249 262 267 304 308 320 321 323 326 363 367 369 370 371 372 373 374 375 376 377 377 378 379 380 381 382 383 384 385 386 387 388 389 390 391 392 393 396 404 408 412 421 431 455 458 459 460 461 468 469 470 471 476 477 478 479 483" ,
	  		template='Vec3d', src ='@topo')	



	

	modelVisu1 = cilinder1.addChild('visu')
	modelVisu1.addObject('MeshSTLLoader', name='loader', filename='mesh/cylinder5296.stl')
	modelVisu1.addObject('OglModel', src='@loader', color=[1,0,0,1])
	modelVisu1.addObject('BarycentricMapping')




## Kelvin-Voigt Cylinder: Material F5000- R0.5  (GREEN)
	cilinder2 = rootNode.addChild('cilinder2')

	cilinder2.addObject('EulerImplicitSolver', name="Solver")

	cilinder2.addObject('CGLinearSolver', name="ItSolver", iterations="250", tolerance="1e-15", threshold = '1e-30')
	cilinder2.addObject('MeshVTKLoader', name='loader', filename='mesh/cylinder1513.vtk', translation = [0, 0.1, 0])
	cilinder2.addObject('MechanicalObject', name='tetras', template='Vec3d', src = '@loader')
	cilinder2.addObject('TetrahedronSetTopologyContainer', name="topo", src ='@loader')
	cilinder2.addObject('TetrahedronSetTopologyModifier' ,  name="Modifier")
	cilinder2.addObject('TetrahedronSetGeometryAlgorithms', template="Vec3d" ,name="GeomAlgo")

	cilinder2.addObject('UniformMass', totalMass="0.0126", src = '@topo')

	cilinder2.addObject('BoxROI', name='boxROI1',box="-0.011 0.09 -0.001  0.011 0.111 0.001", drawBoxes=True)
	cilinder2.addObject('FixedConstraint', indices = '@boxROI1.indices')
	E1 = 1068000
	tau1 = 0.897
	k0 = 1/1.150e-3
	cilinder2.addObject('TetrahedronViscoelasticityFEMForceField', template='Vec3d', name='FEM', src ='@topo',materialName="KelvinVoigtFirstOrder", ParameterSet= str(E1)+' '+str(tau1)+' '+str(k0))

	cilinder2.addObject('ConstantForceField', name = "CFF", listening = True, force =[0,0,0],indices ="4 5 6 7 20 21 22 23 28 29 30 31 68 69 87 88 106 107 143 144 146 149 203 205 249 262 267 304 308 320 321 323 326 363 367 369 370 371 372 373 374 375 376 377 377 378 379 380 381 382 383 384 385 386 387 388 389 390 391 392 393 396 404 408 412 421 431 455 458 459 460 461 468 469 470 471 476 477 478 479 483" ,
	  		template='Vec3d', src ='@topo')	



	modelVisu2 = cilinder2.addChild('visu')
	modelVisu2.addObject('MeshSTLLoader', name='loader', filename='mesh/cylinder5296.stl', translation = [0,0.1,0])
	modelVisu2.addObject('OglModel', src='@loader', color=[0,1,0,1])
	modelVisu2.addObject('BarycentricMapping')



## SLS-MAXWELL FIRST ORDER Cylinder: Material F5000- R0.5  (BLUE)
	cilinder3 = rootNode.addChild('cilinder3')

	cilinder3.addObject('EulerImplicitSolver', name="Solver")

	cilinder3.addObject('CGLinearSolver', name="ItSolver", iterations="250", tolerance="1e-15", threshold = '1e-30')
	cilinder3.addObject('MeshVTKLoader', name='loader', filename='mesh/cylinder1513.vtk', translation = [0, 0.2, 0])
	cilinder3.addObject('MechanicalObject', name='tetras', template='Vec3d', src = '@loader')
	cilinder3.addObject('TetrahedronSetTopologyContainer', name="topo", src ='@loader')
	cilinder3.addObject('TetrahedronSetTopologyModifier' ,  name="Modifier")
	cilinder3.addObject('TetrahedronSetGeometryAlgorithms', template="Vec3d" ,name="GeomAlgo")

	cilinder3.addObject('UniformMass', totalMass="0.0126", src = '@topo')

	cilinder3.addObject('BoxROI', name='boxROI1',box="-0.011 0.19 -0.001  0.011 0.211 0.001", drawBoxes=True)
	cilinder3.addObject('FixedConstraint', indices = '@boxROI1.indices')
	Einf = 725000
	E1 = 508000
	tau1 = 0.897 
	k0 = 1/1.150e-3
	cilinder3.addObject('TetrahedronViscoelasticityFEMForceField', template='Vec3d', name='FEM', src ='@topo',materialName="SLSMaxwellFirstOrder", ParameterSet= str(Einf)+' '+str(E1)+' '+str(tau1)+' '+str(k0))

	cilinder3.addObject('ConstantForceField', name = "CFF", listening = True, force =[0,0,0],indices ="4 5 6 7 20 21 22 23 28 29 30 31 68 69 87 88 106 107 143 144 146 149 203 205 249 262 267 304 308 320 321 323 326 363 367 369 370 371 372 373 374 375 376 377 377 378 379 380 381 382 383 384 385 386 387 388 389 390 391 392 393 396 404 408 412 421 431 455 458 459 460 461 468 469 470 471 476 477 478 479 483" ,
	  		template='Vec3d', src ='@topo')	



	modelVisu3 = cilinder3.addChild('visu')
	modelVisu3.addObject('MeshSTLLoader', name='loader', filename='mesh/cylinder5296.stl', translation = [0,0.2,0])
	modelVisu3.addObject('OglModel', src='@loader', color=[0,0,1,1])
	modelVisu3.addObject('BarycentricMapping')



## SLS-KELVIN-VOIGT FIRST ORDER Cylinder: Material F5000- R0.5  (Cyan)
	cilinder4 = rootNode.addChild('cilinder4')

	cilinder4.addObject('EulerImplicitSolver', name="Solver")

	cilinder4.addObject('CGLinearSolver', name="ItSolver", iterations="250", tolerance="1e-15", threshold = '1e-30')
	cilinder4.addObject('MeshVTKLoader', name='loader', filename='mesh/cylinder1513.vtk', translation = [0, 0.3, 0])
	cilinder4.addObject('MechanicalObject', name='tetras', template='Vec3d', src = '@loader')
	cilinder4.addObject('TetrahedronSetTopologyContainer', name="topo", src ='@loader')
	cilinder4.addObject('TetrahedronSetTopologyModifier' ,  name="Modifier")
	cilinder4.addObject('TetrahedronSetGeometryAlgorithms', template="Vec3d" ,name="GeomAlgo")

	cilinder4.addObject('UniformMass', totalMass="0.0126", src = '@topo')

	cilinder4.addObject('BoxROI', name='boxROI1',box="-0.011 0.29 -0.001  0.011 0.311 0.001", drawBoxes=True)
	cilinder4.addObject('FixedConstraint', indices = '@boxROI1.indices')

	Einf = 725000
	E1 = 508000
	tau1 = 0.897
	k0 = 1/1.150e-3
	cilinder4.addObject('TetrahedronViscoelasticityFEMForceField', template='Vec3d', name='FEM', src ='@topo',materialName="SLSKelvinVoigtFirstOrder", ParameterSet= str(Einf)+' '+str(Einf)+' '+str(tau1)+' '+str(k0))

	cilinder4.addObject('ConstantForceField', name = "CFF", listening = True, force =[0,0,0],indices ="4 5 6 7 20 21 22 23 28 29 30 31 68 69 87 88 106 107 143 144 146 149 203 205 249 262 267 304 308 320 321 323 326 363 367 369 370 371 372 373 374 375 376 377 377 378 379 380 381 382 383 384 385 386 387 388 389 390 391 392 393 396 404 408 412 421 431 455 458 459 460 461 468 469 470 471 476 477 478 479 483" ,
	  		template='Vec3d', src ='@topo')	


	cilinder4.addObject(CylinderController(node=rootNode, CFF1 = rootNode.cilinder1.CFF, CFF2 = rootNode.cilinder2.CFF, CFF3 = rootNode.cilinder3.CFF, CFF4 = rootNode.cilinder4.CFF, pos1 = rootNode.cilinder1.tetras.position.value, pos2 = rootNode.cilinder2.tetras.position.value, pos3 = rootNode.cilinder3.tetras.position.value, pos4 = rootNode.cilinder4.tetras.position.value))

	modelVisu4 = cilinder4.addChild('visu')
	modelVisu4.addObject('MeshSTLLoader', name='loader', filename='mesh/cylinder5296.stl', translation = [0,0.3,0])
	modelVisu4.addObject('OglModel', src='@loader', color=[0,1,1,1])
	modelVisu4.addObject('BarycentricMapping')



	return rootNode
