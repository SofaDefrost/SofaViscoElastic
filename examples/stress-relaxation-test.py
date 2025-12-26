import os
import math
import numpy as np

import Sofa.Core
# import SofaRuntime  # Not required if you use runSofa with SofaPython3 plugin


# Base path of this script (used to build robust relative paths)
BASE_DIR = os.path.dirname(os.path.abspath(__file__))


class CylinderController(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        # Base class init
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        self.time = 0.0
        self.node = kwargs["node"]
        # Positions of the mechanical DOFs at initialization
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

    def onAnimateEndEvent(self, event):
        # Current simulation time
        self.time = self.node.time.value

        # Access Cauchy stress field on FEM (array of Vec3)
        stress_field = self.node.cylinder.FEM.CauchyStress.value
        if len(stress_field) == 0:
            return

        # Stress at max-Z node
        stress = stress_field[self.posmax1]

        # Optional: parse relaxation parameters from ParameterSet if needed
        # ParameterSet is a string like "G1 tau1 G2 tau2 lamb"
        try:
            params = list(map(float, self.node.cylinder.FEM.ParameterSet.value.split()))
            self.tau = params[1] if len(params) > 1 else 0.0
        except Exception:
            self.tau = 0.0

        # Engineering strain at max-Z node
        epsilon = (
            self.node.cylinder.tetras.position.value[self.posmax1][2] - self.lin
        ) / self.lin

        # Print Z-stress in MPa
        print(
            f"t = {self.time:.6g} s, epsilon = {epsilon*100:.4f} %, "
            f"sigma_z = {stress[2] / 1e6:.4f} MPa"
        )


# In this code we do a stress relaxation test:
# we apply a displacement step using a PositionConstraint.


def createScene(rootNode):

    # SofaPython3 does not register components → RequiredPlugin is not strictly needed
    # rootNode.addObject('RequiredPlugin', name='SofaPython3')

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
    rootNode.addObject("RequiredPlugin", name="SofaViscoElastic")

    # Needed for PositionConstraint (SoftRobots plugin)
    rootNode.addObject("RequiredPlugin", name="SoftRobots")

    # Animation / constraint loop
    rootNode.addObject("FreeMotionAnimationLoop")

    # BlockGaussSeidelConstraintSolver was removed in SOFA 25.06 → use GenericConstraintSolver
    rootNode.addObject(
        "GenericConstraintSolver",
        maxIterations=1e4,
        tolerance=1e-12,
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

    # Robust path for the volumetric mesh
    vtk_path = os.path.join(BASE_DIR, "mesh", "cylinder1513.vtk")

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

    # Viscoelastic material parameters for Burgers model
    G1 = 70e6
    tau1 = 1e6 / G1
    G2 = 20e6
    tau2 = 1e6 / G2
    G3 = 10e6
    tau3 = 1e6 / G3
    nu = 0.44
    # lambda = 2 * G1 * nu / (1 - 2 * nu)  # Lame parameter (if needed)
    lamb = 0.0

    cylinder.addObject(
        "TetrahedronViscoelasticityFEMForceField",
        template="Vec3d",
        name="FEM",
        src="@topo",
        materialName="Burgers",
        ParameterSet=f"{G1} {tau1} {G2} {tau2} {lamb}",
    )

    # Fixed boundary region
    cylinder.addObject(
        "BoxROI",
        name="boxROI",
        box="-0.011 -0.011 -0.001  0.011 0.011 0.001",
        drawBoxes=True,
    )
    cylinder.addObject("FixedProjectiveConstraint", indices="@boxROI.indices")

    # Region for imposed displacement
    cylinder.addObject(
        "BoxROI",
        name="boxToPull",
        box=[-0.011, -0.011, 0.1, 0.011, 0.011, 0.101],
        drawBoxes=False,
    )
    cylinder.addObject(
        "PartialFixedProjectiveConstraint",
        indices=cylinder.boxToPull.indices.linkpath,
        fixedDirections=[1, 1, 0],
    )

    # Step input using PositionConstraint (displacement in Z)
    cylinder.addObject(
        "PositionConstraint",
        name="displacement",
        indices=cylinder.boxToPull.indices.linkpath,
        valueType="displacement",
        value=1e-2,
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

    # Visual model
    stl_path = os.path.join(BASE_DIR, "mesh", "cylinder5296.stl")
    modelVisu3 = cylinder.addChild("visu")
    modelVisu3.addObject(
        "MeshSTLLoader",
        name="loader",
        filename=stl_path,
        translation=[0.0, 0.0, 0.0],
    )
    modelVisu3.addObject("OglModel", src="@loader", color=[1.0, 1.0, 1.0, 1.0])
    modelVisu3.addObject("BarycentricMapping")

    return rootNode
