# to be able to add sofa objects you need to first load the plugins that implement them.
# For simplicity you can load the plugin "SofaComponentAll" that will load all most
# common sofa objects.
import SofaRuntime
# SofaRuntime.importPlugin("SofaComponentAll")

# to add elements like Node or objects
import Sofa.Core
root = Sofa.Core.Node()

import math 
import numpy as np
import matplotlib.pyplot as plt

import os
path = os.path.dirname(os.path.abspath(__file__))+'/plot/'


class PendulumController(Sofa.Core.Controller):

	def __init__(self, *args, **kwargs):
		Sofa.Core.Controller.__init__(self,*args, **kwargs) #needed
		self.node = kwargs['node']		
		# self.Positions =  kwargs['position']
		# file3 = open("/data/Softwares/sofa/src/master/forum/Pasquale/plot/EulerExplicit.txt","w")
		# file3.write(str(0.0)+' '+str(self.Positions[1][0])+' '+str(self.Positions[1][1])+' '+ str(self.Positions[1][2])+ ' '+'\n' )
		# file3.close()
		self.times = []
		self.resultEE = []
		self.resultRK = []
		self.resultHHT = []
		self.resultNewmark = []
		self.resultEI = []
		self.resultTRAP = []


		file1 = open("/home/pasquale/EulerExplicit.txt","w")
		file1.write(str(0.0)+' '+str(self.node.EulerExplicit.Particles.position.value[1][2])+'\n')
		file1.close()

		file2 = open("/home/pasquale/RungeKutta.txt","w")
		file2.write(str(0.0)+' '+str(self.node.RungeKutta.Particles.position.value[1][2])+'\n')
		file2.close()


		file3 = open("/home/pasquale/HHT.txt","w")
		file3.write(str(0.0)+' '+str(self.node.HHT.Particles.position.value[1][2])+'\n')
		file3.close()

		file4 = open("/home/pasquale/Newmark.txt","w")
		file4.write(str(0.0)+' '+str(self.node.Newmark.Particles.position.value[1][2])+'\n')
		file4.close()

		file5 = open("/home/pasquale/EulerImplicit.txt","w")
		file5.write(str(0.0)+' '+str(self.node.EulerImplicit.Particles.position.value[1][2])+'\n')
		file5.close()

		file5 = open("/home/pasquale/Trapezoidal.txt","w")
		file5.write(str(0.0)+' '+str(self.node.Trapezoidal.Particles.position.value[1][2])+'\n')
		file5.close()


	def onAnimateBeginEvent(self, event): # called at each begin of animation step

		self.time = self.node.time.value

		file1 = open("/home/pasquale/EulerExplicit.txt","a")
		file1.write(str(self.time)+' '+str(self.node.EulerExplicit.Particles.position.value[1][2])+'\n')
		file1.close()

		file2 = open("/home/pasquale/RungeKutta.txt","a")
		file2.write(str(self.time)+' '+str(self.node.RungeKutta.Particles.position.value[1][2])+'\n')
		file2.close()


		file3 = open("/home/pasquale/HHT.txt","a")
		file3.write(str(self.time)+' '+str(self.node.HHT.Particles.position.value[1][2])+'\n')
		file3.close()

		file4 = open("/home/pasquale/Newmark.txt","a")
		file4.write(str(self.time)+' '+str(self.node.Newmark.Particles.position.value[1][2])+'\n')
		file4.close()

		file5 = open("/home/pasquale/EulerImplicit.txt","a")
		file5.write(str(self.time)+' '+str(self.node.EulerImplicit.Particles.position.value[1][2])+'\n')
		file5.close()

		file5 = open("/home/pasquale/Trapezoidal.txt","a")
		file5.write(str(self.time)+' '+str(self.node.Trapezoidal.Particles.position.value[1][2])+'\n')
		file5.close()

def createScene(rootNode):


	rootNode.addObject('RequiredPlugin', name="Sofa.Component.Visual")
	rootNode.addObject('RequiredPlugin', name="Sofa.Component.LinearSolver.Iterative")
	rootNode.addObject('RequiredPlugin', name="Sofa.Component.ODESolver.Backward")
	rootNode.addObject('RequiredPlugin', name="Sofa.Component.Setting")
	rootNode.addObject('RequiredPlugin', name="Sofa.GL.Component.Rendering3D")
	rootNode.addObject('RequiredPlugin', name="Sofa.Component.Constraint.Projective")
	rootNode.addObject('RequiredPlugin', name="Sofa.Component.Mass")
	rootNode.addObject('RequiredPlugin', name="Sofa.Component.SolidMechanics.Spring")
	rootNode.addObject('RequiredPlugin', name="Sofa.Component.StateContainer")
	rootNode.addObject('RequiredPlugin', name="MyAwesomeComponents")
	rootNode.addObject('RequiredPlugin', name="Sofa.Component.Collision.Geometry")
	rootNode.addObject('RequiredPlugin', name="Sofa.Component.ODESolver.Forward")
	rootNode.addObject('VisualStyle', displayFlags='hideVisualModels showBehaviorModels showCollisionModels hideBoundingCollisionModels hideForceFields hideInteractionForceFields hideWireframe')

	rootNode.gravity=[ 0, 0, 9.81]
	rootNode.dt = 0.01
	rootNode.name = 'rootNode'
	rootNode.addObject('DefaultAnimationLoop', computeBoundingBox="0")

	rootNode.addObject('BackgroundSetting', color=[0, 0.168627, 0.211765, 1])
	rootNode.addObject('OglSceneFrame', style='Arrows', alignment='TopRight')

	EE = rootNode.addChild("EulerExplicit")
	EE.addObject('EulerExplicitSolver', name="Solver")
	EE.addObject('CGLinearSolver', name="CG Solver", iterations="25", tolerance="1e-10", threshold="1e-10")
	EE.addObject('MechanicalObject', name="Particles", template="Vec3d",position="0 0 0 1 0 0", velocity="0 0 0 0 0 0", showObject = True, showObjectScale = 10)
	EE.addObject('UniformMass', name="Mass", totalMass="0.1")
	EE.addObject('FixedConstraint', indices="0")
	EE.addObject('StiffSpringForceField', name="Springs", spring="0 1 100 0.0 1")
	# EE.addObject('SphereCollisionModel', radius="0.1")



	RK = rootNode.addChild("RungeKutta")
	RK.addObject('RungeKutta4Solver', name="Solver")
	RK.addObject('CGLinearSolver', name="CG Solver", iterations="25", tolerance="1e-10", threshold="1e-10")
	RK.addObject('MechanicalObject', name="Particles", template="Vec3d",position="0 0 0 1 0 0", velocity="0 0 0 0 0 0", showObject = True, showObjectScale = 10)
	RK.addObject('UniformMass', name="Mass", totalMass="0.1")
	RK.addObject('FixedConstraint', indices="0")
	RK.addObject('StiffSpringForceField', name="Springs", spring="0 1 100 0.0 1")
	# RK.addObject('SphereCollisionModel', radius="0.1")


	RK = rootNode.addChild("HHT")
	RK.addObject('HHTSolver', name="Solver")
	RK.addObject('CGLinearSolver', name="CG Solver", iterations="25", tolerance="1e-10", threshold="1e-10")
	RK.addObject('MechanicalObject', name="Particles", template="Vec3d",position="0 0 0 1 0 0", velocity="0 0 0 0 0 0", showObject = True, showObjectScale = 10)
	RK.addObject('UniformMass', name="Mass", totalMass="0.1")
	RK.addObject('FixedConstraint', indices="0")
	RK.addObject('StiffSpringForceField', name="Springs", spring="0 1 100 0.0 1")
	# RK.addObject('SphereCollisionModel', radius="0.1")



	BDF = rootNode.addChild("Newmark")
	BDF.addObject('NewmarkImplicitSolver', name="Solver")
	BDF.addObject('CGLinearSolver', name="CG Solver", iterations="25", tolerance="1e-10", threshold="1e-10")
	BDF.addObject('MechanicalObject', name="Particles", template="Vec3d",position="0 0 0 1 0 0", velocity="0 0 0 0 0 0", showObject = True, showObjectScale = 10)
	BDF.addObject('UniformMass', name="Mass", totalMass="0.1")
	BDF.addObject('FixedConstraint', indices="0")
	BDF.addObject('StiffSpringForceField', name="Springs", spring="0 1 100 0.0 1")
	# BDF.addObject('SphereCollisionModel', radius="0.1")



	EI = rootNode.addChild("EulerImplicit")
	EI.addObject('EulerImplicitSolver', name="Solver")
	EI.addObject('CGLinearSolver', name="CG Solver", iterations="25", tolerance="1e-10", threshold="1e-10")
	EI.addObject('MechanicalObject', name="Particles", template="Vec3d",position="0 0 0 1 0 0", velocity="0 0 0 0 0 0", showObject = True, showObjectScale = 10)
	EI.addObject('UniformMass', name="Mass", totalMass="0.1")
	EI.addObject('FixedConstraint', indices="0")
	EI.addObject('StiffSpringForceField', name="Springs", spring="0 1 100 0.0 1")
	# EI.addObject('SphereCollisionModel', radius="0.1")


	TRAP = rootNode.addChild("Trapezoidal")
	TRAP.addObject('EulerImplicitSolver', name="Solver", trapezoidalScheme=True)
	TRAP.addObject('CGLinearSolver', name="CG Solver", iterations="25", tolerance="1e-10", threshold="1e-10")
	TRAP.addObject('MechanicalObject', name="Particles", template="Vec3d",position="0 0 0 1 0 0", velocity="0 0 0 0 0 0", showObject = True, showObjectScale = 10)
	TRAP.addObject('UniformMass', name="Mass", totalMass="0.1")
	TRAP.addObject('FixedConstraint', indices="0")
	TRAP.addObject('StiffSpringForceField', name="Springs", spring="0 1 100 0.0 1")
	# TRAP.addObject('SphereCollisionModel', radius="0.1")


	rootNode.addObject(PendulumController(node=rootNode))

	return rootNode
