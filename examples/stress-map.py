import os
import math
import numpy as np

import Sofa.Core
# import SofaRuntime  # Not strictly needed if you launch with runSofa + SofaPython3 plugin


# Base path of this script (used for relative mesh paths)
BASE_DIR = os.path.dirname(os.path.abspath(__file__))


class CylinderController(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        self.time = 0.0
        self.node = kwargs["node"]
        # Expect positions as a numpy array / list of Vec3
        self.pos = kwargs["pos"]

        self.max1 = 0.0
        self.posmax1 = 0

        # Find node with maximum Z coordinate
        for j in range(len(self.pos)):
            if self.pos[j][2] >= self.max1:
                self.max1 = self.pos[j][2]
                self.posmax1 = j

        print("Max-Z node: z =", self.pos[self.posmax1][2], "index =", self.posmax1)

        # Reference length for strain computation
        self.lin = self.node.cylinder.tetras.position.value[self.posmax1][2]

    def onAnimateBeginEvent(self, event):
        # CauchyStress is assumed in Pa → divide by 1e6 to get MPa
        self.stress = self.node.cylinder.FEM.CauchyStress.value / 1e6  # MPa

        if len(self.stress) == 0:
            return

        self.time = self.node.time.value

        # Strain at the max-Z node
        epsilon = (
            self.node.cylinder.tetras.position.value[self.posmax1][2] - self.lin
        ) / self.lin

        # Map stress (Z-component) to visualization data
        self.node.cylinder.visu.display.pointData.value = self.stress[0 : len(self.pos), 2]
        self.node.cylinder.visu.Map.min.value = (
            self.node.cylinder.visu.display.currentMin.value
        )
        self.node.cylinder.visu.Map.max.value = (
            self.node.cylinder.visu.display.currentMax.value
        )

        # self.stress is already in MPa
        print("Stress at max-Z node [MPa] =", self.stress[self.posmax1][2])


# In this code we perform a stress relaxation test:
# we apply a displacement step using a PositionConstraint.


def createScene(rootNode):
    # SofaPython3 bindings do not register BaseObject → RequiredPlugin is not needed
    # rootNode.addObject("RequiredPlugin", name="SofaPython3")

    rootNode.addObject("RequiredPlugin", name="Sofa.Component.Engine.Select")
    rootNode.addObject("RequiredPlugin", name="Sofa.Component.IO.Mesh")
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
    rootNode.addObject("RequiredPlugin", name="Sofa.GL.Component.Rendering2D")  # For OglColorMap
    rootNode.addObject("RequiredPlugin", name="SoftRobots")  # For PositionConstraint
    rootNode.addObject("RequiredPlugin", name="SofaViscoElastic")

    # Animation / constraints loop
    rootNode.addObject("FreeMotionAnimationLoop")

    # BlockGaussSeidelConstraintSolver was removed in SOFA 25.06 → use GenericConstraintSolver
    rootNode.addObject(
        "GenericConstraintSolver",
        maxIterations=1e4,
        tolerance=1e-50,
    )

    rootNode.gravity = [0.0, 0.0, -9.81]
    rootNode.dt = 1e6 / (20e6 * 100)

    rootNode.addObject("VisualStyle", displayFlags="hideForceFields")
    rootNode.addObject("OglSceneFrame", style="Arrows", alignment="TopRight")

    # ======================
    # Cylinder definition
    # ======================
    cylinder = rootNode.addChild("cylinder")

    cylinder.addObject(
        "EulerImplicitSolver",
        name="Solver",
        rayleighMass=0.0,
        rayleighStiffness=0.0,
        firstOrder=True,
        trapezoidalScheme=False,
    )
    cylinder.addObject("SparseLDLSolver", name="directsolver")

    # Robust path for the mesh file
    vtk_path = os.path.join(BASE_DIR, "mesh", "cylinder9178.vtk")

    cylinder.addObject("MeshVTKLoader", name="loader", filename=vtk_path)
    cylinder.addObject(
        "MechanicalObject",
        name="tetras",
        template="Vec3d",
        src="@loader",
        rotation=[0.0, 0.0, 0.0],
    )
    cylinder.addObject("TetrahedronSetTopologyContainer", name="topo", src="@loader")
    cylinder.addObject("TetrahedronSetTopologyModifier", name="Modifier")
    cylinder.addObject(
        "TetrahedronSetGeometryAlgorithms",
        template="Vec3d",
        name="GeomAlgo",
    )

    cylinder.addObject("UniformMass", totalMass="0.0126", src="@topo")

    # Viscoelastic material parameters
    G1 = 70e6
    tau1 = 1e6 / G1
    G2 = 20e6
    tau2 = 1e6 / G2
    G3 = 10e6
    tau3 = 1e9 / G3
    nu = 0.44
    lamb = 30e6

    cylinder.addObject(
        "TetrahedronViscoelasticityFEMForceField",
        template="Vec3d",
        name="FEM",
        src="@topo",
        materialName="MaxwellFirstOrder",
        ParameterSet=f"{G1} {tau1} {lamb}",
    )

    # Fixed boundary region
    cylinder.addObject(
        "BoxROI",
        name="boxROI",
        box="-0.011 -0.011 -0.001  0.011 0.011 0.001",
        drawBoxes=True,
    )
    cylinder.addObject("FixedProjectiveConstraint", indices="@boxROI.indices")

    # ROI for imposed displacement
    cylinder.addObject(
        "BoxROI",
        name="boxToPull",
        box=[-0.011, -0.011, 0.1, 0.011, 0.011, 0.101],
        drawBoxes=True,
    )
    cylinder.addObject(
        "PartialFixedProjectiveConstraint",
        indices=cylinder.boxToPull.indices.linkpath,
        fixedDirections=[1, 1, 0],
    )

    # Step input using PositionConstraint
    cylinder.addObject(
        "PositionConstraint",
        name="displacement",
        indices=cylinder.boxToPull.indices.linkpath,
        valueType="displacement",
        value=1e-4,
        useDirections=[0, 0, 1],
    )
    cylinder.addObject("LinearSolverConstraintCorrection")

    # Attach Python controller
    cylinder.addObject(
        CylinderController(
            node=rootNode,
            pos=rootNode.cylinder.tetras.position.value,
        )
    )

    # Visualization: stress mapped through DataDisplay + OglColorMap
    modelVisu1 = cylinder.addChild("visu")
    modelVisu1.addObject("DataDisplay", name="display")
    modelVisu1.addObject(
        "OglColorMap",
        name="Map",
        colorScheme="Blue to Red",
        showLegend=True,
    )
    modelVisu1.addObject("IdentityMapping", input="@..", output="@.")

    return rootNode
