# Python scene for SOFA 25.06 – creep / amplitude sweep test on visco-hyperelastic cylinder

import os
import math
import numpy as np
# from scipy import signal   # Enable only if SciPy is installed

import Sofa.Core

# Base path of this script (used for mesh paths)
BASE_DIR = os.path.dirname(os.path.abspath(__file__))


class CylinderController(Sofa.Core.Controller):

    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.time = 0.0
        self.node = kwargs['node']
        self.pos3 = kwargs['pos3']

        self.max1 = 0.0
        self.posmax1 = 0

        # Find node with maximum Z coordinate
        for j in range(len(self.pos3)):
            if self.pos3[j][2] >= self.max1:
                self.max1 = self.pos3[j][2]
                self.posmax1 = j

        print("Index of max-Z node:", self.posmax1)

        # Reference length for strain computation
        self.lin = self.node.cylinder.tetras.position.value[self.posmax1][2]

    def onAnimateBeginEvent(self, event):

        # Current simulation time
        self.time = self.node.time.value

        # Safe extraction of relaxation time (if needed)
        try:
            self.tau = float(self.node.cylinder.FEM.ParameterSet.value.split()[2])
        except Exception:
            self.tau = 0.0

        pos = self.node.cylinder.tetras.position.value

        # Strain computed on node index 4 (check for safety)
        if len(pos) > 4:
            epsilon = (pos[4][2] - self.lin) / self.lin
        else:
            epsilon = 0.0

        # ===== LOAD INPUT =====
        # Step input in Z direction
        self.node.cylinder.CFF.totalForce.value = [
            0.0,
            0.0,
            (3.141592653589793e2 / 2.0) * np.heaviside(self.time, 0.0),
        ]

        # You can extend this to amplitude sweep using time intervals and tau

        print("epsilon [%] =", epsilon * 100.0)


def createScene(rootNode):

    rootNode.addObject(
        'VisualStyle',
        displayFlags='showVisualModels hideBehaviorModels hideCollisionModels '
                     'hideBoundingCollisionModels hideForceFields '
                     'hideInteractionForceFields hideWireframe'
    )

    # SofaPython3 bindings do not register BaseObject → RequiredPlugin is not needed
    # rootNode.addObject('RequiredPlugin', name='SofaPython3')

    # Required core Sofa component plugins
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
    rootNode.addObject("RequiredPlugin", name="Sofa.Component.Constraint.Projective")
    rootNode.addObject("RequiredPlugin", name="Sofa.Component.ODESolver.Backward")
    rootNode.addObject("RequiredPlugin", name="SofaViscoElastic")

    # Global simulation settings
    rootNode.gravity = [0, 0, -9.81]
    rootNode.dt = 1e6 / (10e6 * 100)
    rootNode.name = 'rootNode'
    rootNode.addObject('DefaultAnimationLoop', computeBoundingBox="0")

    # BlockGaussSeidelConstraintSolver was removed in 25.06 → use GenericConstraintSolver
    rootNode.addObject('GenericConstraintSolver', tolerance=1e-24, maxIterations=100000000)

    rootNode.addObject('OglSceneFrame', style='Arrows', alignment='TopRight')

    # =========================
    #  VISCO-HYPERELASTIC CYLINDER
    # =========================
    cylinder = rootNode.addChild('cylinder')

    cylinder.addObject(
        'EulerImplicitSolver',
        name="Solver",
        rayleighMass=0.0,
        rayleighStiffness=0.0
    )

    cylinder.addObject(
        'CGLinearSolver',
        name="ItSolver",
        iterations="250",
        tolerance="1e-12",
        threshold='1e-12'
    )

    # Mesh files (relative path)
    vtk_path = os.path.join(BASE_DIR, 'mesh', 'cylinder5296.vtk')
    stl_path = os.path.join(BASE_DIR, 'mesh', 'cylinder5296.stl')

    cylinder.addObject('MeshVTKLoader', name='loader',
                       filename=vtk_path, translation=[0, 0.0, 0])

    cylinder.addObject('MechanicalObject', name='tetras',
                       template='Vec3d', src='@loader')

    cylinder.addObject('TetrahedronSetTopologyContainer',
                       name="topo", src='@loader')

    cylinder.addObject('TetrahedronSetTopologyModifier', name="Modifier")

    cylinder.addObject('TetrahedronSetGeometryAlgorithms',
                       template="Vec3d", name="GeomAlgo")

    cylinder.addObject('UniformMass', totalMass="1", src='@topo')

    # Fixed boundary ROI
    cylinder.addObject('BoxROI', name='boxROI1',
                       box="-0.011 -0.011 -0.001  0.011 0.011 0.001",
                       drawBoxes=True)

    cylinder.addObject('FixedProjectiveConstraint',
                       indices='@boxROI1.indices')

    # Loading ROI
    cylinder.addObject('BoxROI', name="boxToPull",
                       box=[-0.011, -0.011, 0.1, 0.011, 0.011, 0.101],
                       drawBoxes=False)

    # Material parameters
    mu = 70e6
    E2 = 20e6
    E3 = 10e6
    tau1 = 1e6 / mu
    k = 44e6

    # Visco-hyperelastic force field
    cylinder.addObject(
        'TetrahedronViscoHyperelasticityFEMForceField',
        template='Vec3d',
        name='FEM',
        src='@topo',
        materialName="SLSNeoHookeanFirstOrder",
        ParameterSet=f"{mu} {mu} {tau1} {k}"
    )

    # External force applied on ROI nodes
    cylinder.addObject(
        'ConstantForceField',
        name="CFF",
        listening=True,
        totalForce=[0, 0, 0],
        template="Vec3d",
        src="@topo",
        indices=cylinder.boxToPull.indices.linkpath
    )

    # Python controller
    cylinder.addObject(
        CylinderController(
            node=rootNode,
            pos3=rootNode.cylinder.tetras.position.value
        )
    )

    # Visual model
    modelVisu3 = cylinder.addChild('visu')
    modelVisu3.addObject('MeshSTLLoader', name='loader',
                         filename=stl_path, translation=[0.0, 0.0, 0.0])
    modelVisu3.addObject('OglModel', src='@loader', color=[1, 0, 0, 1])
    modelVisu3.addObject('BarycentricMapping')

    return rootNode
