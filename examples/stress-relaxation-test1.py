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
		self.pos = kwargs['pos']


		self.max1 = 0
		self.posmax1 = 0

		#for j in range(0, len(self.pos)):
		#	if self.pos[j][2] >= self.max1 and self.pos[j][0]>=-0.002 and self.pos[j][0]<= 0.002 and self.pos[j][1]>=-0.002 and self.pos[j][1]<= 0.002:
		#		self.max1 = self.pos[j][2]
		#		self.posmax1 = j
		#print(self.posmax1)

		#file1 = open("/home/pasquale/Script_Sofa/example/plot/posMaxwell-stress-relaxation.txt","w")
		#file1.write(str(0.0)+' '+str(self.pos[369][2])+ '\n' )
		#file1.close()



	def onAnimateBeginEvent(self,event):
		self.time += self.node.dt.value 
		#file1 = open("/home/pasquale/Script_Sofa/example/plot/posMaxwell-stress-relaxation.txt","a")
		#file1.write(str(self.time)+' '+str(self.node.cilinder1.tetras.position.value[369][2])+'\n' )
		#file1.close()


def createScene(rootnode):
	rootnode.addObject('FreeMotionAnimationLoop')
	rootnode.addObject('GenericConstraintSolver', maxIterations=1e4, tolerance=1e-50)
	rootnode.gravity = [0,0,0]
	rootnode.addObject('VisualStyle', displayFlags='hideForceFields')
	rootnode.addObject('OglSceneFrame', style='Arrows', alignment='TopRight')	

	cilinder1 = rootnode.addChild('cilinder1')

	cilinder1.addObject('EulerImplicitSolver', name="Solver")
	cilinder1.addObject('SparseLDLSolver', name="directsolver")
	#cilinder1.addObject('CGLinearSolver', name="ItSolver", iterations="250", tolerance="1e-15", threshold = '1e-30')

	cilinder1.addObject('MeshVTKLoader', name='loader', filename='mesh/cylinder1513.vtk',)
	cilinder1.addObject('MechanicalObject', name='tetras', template='Vec3d', src = '@loader', rotation = [0, 0 , 0])
	cilinder1.addObject('TetrahedronSetTopologyContainer', name="topo", src ='@loader')
	cilinder1.addObject('TetrahedronSetTopologyModifier' ,  name="Modifier")
	cilinder1.addObject('TetrahedronSetGeometryAlgorithms', template="Vec3d" ,name="GeomAlgo")

	cilinder1.addObject('UniformMass', totalMass="0.0126", src = '@topo')
	E1 = 1e8
	tau1 = 1
	nu = 0.44
	cilinder1.addObject('TetrahedronViscoelasticityFEMForceField', template='Vec3d', name='FEM', src ='@topo',materialName="MaxwellFirstOrder", ParameterSet= str(E1)+' '+str(tau1)+' '+str(nu))
	cilinder1.addObject('BoxROI', name='boxROI',box="-0.011 -0.011 -0.001  0.011 0.011 0.001", drawBoxes=True)
	cilinder1.addObject('FixedConstraint', indices = '@boxROI.indices')
	cilinder1.addObject('BoxROI', name="boxToPull", box=[-0.011, -0.011, 0.09, 0.011, 0.011, 0.11], drawBoxes=False)
	cilinder1.addObject('PartialFixedConstraint', indices=cilinder1.boxToPull.indices.linkpath, fixedDirections=[1, 1, 0])
	cilinder1.addObject('PositionConstraint', indices=cilinder1.boxToPull.indices.linkpath,valueType="displacement", value=0.0, useDirections=[0, 0, 1])
	cilinder1.addObject('LinearSolverConstraintCorrection')

	cilinder1.addObject(CylinderController(node=rootnode, pos = rootnode.cilinder1.tetras.position))

	

	modelVisu1 = cilinder1.addChild('visu')
	modelVisu1.addObject('MeshSTLLoader', name='loader', filename='mesh/cylinder5296.stl')
	modelVisu1.addObject('OglModel', src='@loader', color=[1,1,1,1])
	modelVisu1.addObject('BarycentricMapping')
	return rootnode