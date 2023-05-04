import SofaRuntime
SofaRuntime.importPlugin("SofaComponentAll")

# to add elements like Node or objects
import Sofa.Core
root = Sofa.Core.Node()

import math 
import numpy as np

import os
path = os.path.dirname(os.path.abspath(__file__))+'/plot/'

def createScene(rootnode):
	rootnode.addObject('FreeMotionAnimationLoop')
	rootnode.addObject('GenericConstraintSolver', maxIterations=1e4, tolerance=1e-5, printLog=True)
	rootnode.gravity = [0, -9.81, 0]
	rootnode.addObject('VisualStyle', displayFlags='hideForceFields')
	rootnode.addObject('OglSceneFrame', style='Arrows', alignment='TopRight')


	cilinder2 = rootnode.addChild('cilinder2')

	cilinder2.addObject('EulerImplicitSolver', name="Solver")
	cilinder2.addObject('SparseLDLSolver', name="directsolver")

	cilinder2.addObject('MeshVTKLoader', name='loader', filename='mesh/cylinder1513.vtk',translation = [0,0.1,0])
	cilinder2.addObject('MechanicalObject', name='tetras', template='Vec3d', src = '@loader', rotation = [0, 0 , 0])
	cilinder2.addObject('TetrahedronSetTopologyContainer', name="topo", src ='@loader')
	cilinder2.addObject('TetrahedronSetTopologyModifier' ,  name="Modifier")
	cilinder2.addObject('TetrahedronSetGeometryAlgorithms', template="Vec3d" ,name="GeomAlgo")

	cilinder2.addObject('UniformMass', totalMass="0.0126", src = '@topo')
	E1 = 1068000
	tau1 = 0.897
	k0 = 1/1.150e-3
	cilinder2.addObject('TetrahedronViscoelasticityFEMForceField', template='Vec3d', name='FEM', src ='@topo',materialName="KelvinVoigtFirstOrder", ParameterSet= str(E1)+' '+str(tau1)+' '+str(k0))
	cilinder2.addObject('GenericConstraintCorrection')
	cilinder2.addObject('BoxROI', name='boxROI',box="-0.011 0.09 -0.001  0.011 0.111 0.001", drawBoxes=True)
	cilinder2.addObject('FixedConstraint', indices = '@boxROI.indices')
	cilinder2.addObject('BoxROI', name="boxToPull", box=[-0.011, 0.09, 0.09, 0.011, 0.111, 0.11], drawBoxes=False)
	cilinder2.addObject('PartialFixedConstraint', indices=cilinder2.boxToPull.indices.linkpath, fixedDirections=[1, 1, 0])
	cilinder2.addObject('PositionConstraint', indices=cilinder2.boxToPull.indices.linkpath,valueType="displacement", value=0.03, useDirections=[0, 0, 1])


	

	modelVisu2 = cilinder2.addChild('visu')
	modelVisu2.addObject('MeshSTLLoader', name='loader', filename='mesh/cylinder5296.stl', translation = [0,0.1,0])
	modelVisu2.addObject('OglModel', src='@loader', color=[0,1,0,1])
	modelVisu2.addObject('BarycentricMapping')

	cilinder3 = rootnode.addChild('cilinder3')

	cilinder3.addObject('EulerImplicitSolver', name="Solver")
	cilinder3.addObject('SparseLDLSolver', name="directsolver")

	cilinder3.addObject('MeshVTKLoader', name='loader', filename='mesh/cylinder1513.vtk',translation = [0,0.2,0])
	cilinder3.addObject('MechanicalObject', name='tetras', template='Vec3d', src = '@loader', rotation = [0, 0 , 0])
	cilinder3.addObject('TetrahedronSetTopologyContainer', name="topo", src ='@loader')
	cilinder3.addObject('TetrahedronSetTopologyModifier' ,  name="Modifier")
	cilinder3.addObject('TetrahedronSetGeometryAlgorithms', template="Vec3d" ,name="GeomAlgo")

	cilinder3.addObject('UniformMass', totalMass="0.0126", src = '@topo')
	Einf = 725000
	E1 = 508000
	tau1 = 0.897 
	k0 = 1/1.150e-3
	cilinder3.addObject('TetrahedronViscoelasticityFEMForceField', template='Vec3d', name='FEM', src ='@topo',materialName="SLSMaxwellFirstOrder", ParameterSet= str(Einf)+' '+str(E1)+' '+str(tau1)+' '+str(k0))
	cilinder3.addObject('GenericConstraintCorrection')
	cilinder3.addObject('BoxROI', name='boxROI',box="-0.011 0.19 -0.001  0.011 0.211 0.001", drawBoxes=True)
	cilinder3.addObject('FixedConstraint', indices = '@boxROI.indices')
	cilinder3.addObject('BoxROI', name="boxToPull", box=[-0.011, 0.19, 0.09, 0.011, 0.211, 0.11], drawBoxes=False)
	cilinder3.addObject('PartialFixedConstraint', indices=cilinder3.boxToPull.indices.linkpath, fixedDirections=[1, 1, 0])
	cilinder3.addObject('PositionConstraint', indices=cilinder3.boxToPull.indices.linkpath,valueType="displacement", value=0.03, useDirections=[0, 0, 1])


	

	modelVisu3 = cilinder3.addChild('visu')
	modelVisu3.addObject('MeshSTLLoader', name='loader', filename='mesh/cylinder5296.stl', translation = [0,0.2,0])
	modelVisu3.addObject('OglModel', src='@loader', color=[0,0,1,1])
	modelVisu3.addObject('BarycentricMapping')


	cilinder4 = rootnode.addChild('cilinder4')

	cilinder4.addObject('EulerImplicitSolver', name="Solver")
	cilinder4.addObject('SparseLDLSolver', name="directsolver")

	cilinder4.addObject('MeshVTKLoader', name='loader', filename='mesh/cylinder1513.vtk',translation = [0,0.3,0])
	cilinder4.addObject('MechanicalObject', name='tetras', template='Vec3d', src = '@loader', rotation = [0, 0 , 0])
	cilinder4.addObject('TetrahedronSetTopologyContainer', name="topo", src ='@loader')
	cilinder4.addObject('TetrahedronSetTopologyModifier' ,  name="Modifier")
	cilinder4.addObject('TetrahedronSetGeometryAlgorithms', template="Vec3d" ,name="GeomAlgo")

	cilinder4.addObject('UniformMass', totalMass="0.0126", src = '@topo')
	Einf = 725000
	E1 = 508000
	tau1 = 0.897 
	k0 = 1/1.150e-3
	cilinder4.addObject('TetrahedronViscoelasticityFEMForceField', template='Vec3d', name='FEM', src ='@topo',materialName="SLSKelvinVoigtFirstOrder", ParameterSet= str(Einf)+' '+str(E1)+' '+str(tau1)+' '+str(k0))
	cilinder4.addObject('GenericConstraintCorrection')
	cilinder4.addObject('BoxROI', name='boxROI',box="-0.011 0.29 -0.001  0.011 0.311 0.001", drawBoxes=True)
	cilinder4.addObject('FixedConstraint', indices = '@boxROI.indices')
	cilinder4.addObject('BoxROI', name="boxToPull", box=[-0.011, 0.29, 0.09, 0.011, 0.311, 0.11], drawBoxes=False)
	cilinder4.addObject('PartialFixedConstraint', indices=cilinder4.boxToPull.indices.linkpath, fixedDirections=[1, 1, 0])
	cilinder4.addObject('PositionConstraint', indices=cilinder4.boxToPull.indices.linkpath,valueType="displacement", value=0.03, useDirections=[0, 0, 1])


	

	modelVisu4 = cilinder4.addChild('visu')
	modelVisu4.addObject('MeshSTLLoader', name='loader', filename='mesh/cylinder5296.stl', translation = [0,0.3,0])
	modelVisu4.addObject('OglModel', src='@loader', color=[0,1,1,1])
	modelVisu4.addObject('BarycentricMapping')
	return rootnode