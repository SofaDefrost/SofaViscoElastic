import os
import math
import numpy as np

import Sofa.Core
# import SofaRuntime  # Not required when using runSofa with SofaPython3


# Base path of this script (used to build robust relative paths)
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Folder where time–strain–stress data will be written
PLOT_DIR = os.path.join(BASE_DIR, "plot")
if not os.path.isdir(PLOT_DIR):
    if os.path.isfile(PLOT_DIR):
        raise ValueError(f"path {PLOT_DIR} already exists and is a file")
    else:
        os.mkdir(PLOT_DIR)


class CylinderController(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        # Base class init
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        self.time = 0.0
        self.node = kwargs["node"]
        # Initial positions
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

        # Initialize output file
        out_path = os.path.join(PLOT_DIR, "SLS_Maxwell_relaxation.txt")
        with open(out_path, "w") as f:
            # time [s], strain [%], stress [MPa]
            f.write(f"{0.0} {0.0} {0.0}\n")

    def onAnimateEndEvent(self, event):
        # Current simulation time
        self.time = self.node.time.value

        # Access second Piola–Kirchhoff stress (stressSPK) at max-Z node
        stress_field = self.node.cylinder.FEM.stressSPK.value
        if len(stress_field) == 0:
            return

        stress = stress_field[self.posmax1]

        # Optional: parse relaxation parameters from ParameterSet if needed
        try:
            params = list(map(float, self.node.cylinder.FEM.ParameterSet.value.split()))
            # Here tau is taken from the 3rd entry (index 2) if present
            self.tau = params[2] if len(params) > 2 else 0.0
        except Exception:
            self.tau = 0.0

        # Engineering strain at max-Z node
        epsilon = (
            self.node.cylinder.tetras.position.value[self.posmax1][2] - self.lin
        ) / self.lin

        # Print stress in MPa
        print(
            f"t = {self.time:.6g} s, epsilon = {epsilon*100:.4f} %, "
            f"SPK_z = {stress[2] / 1e6:.4f} MPa"
        )

        # In this code we perform a stress relaxation test:
        # we apply a displacement step using a PositionConstraint.
        out_path = os.path.join(PLOT_DIR, "SLS_Maxwell_relaxation.txt")
        with open(out_path, "a") as f:
            # time [s], strain [%], stress [MPa]
            f.write(f"{self.time} {epsilon*100} {stress[2]/1e6}\n")


def createScene(rootNode):

    # SofaPython3 does not register BaseObject → RequiredPlugin is not strictly needed
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
    rootNode.addObject("RequiredPlugin", name="SoftRobots")       # For PositionConstraint
    rootNode.addObject("RequiredPlugin", name="SofaViscoElastic") # For visco-hyperelastic components

    # Animation / constraint loop
    rootNode.addObject("FreeMotionAnimationLoop")

    # BlockGaussSeidelConstraintSolver was removed in SOFA 25.06 → use GenericConstraintSolver
    rootNode.addObject(
        "GenericConstraintSolver",
        maxIterations=1e4,
        tolerance=1e-50,
    )

    # Gravity in mm units (if you use mm) → -9810 mm/s² along Y
    rootNode.gravity = [0.0, -9810.0, 0.0]
    rootNode.dt = 1e9 / (20e9 * 100)

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

    # Visco-hyperelastic material parameters (SLS Neo-Hookean 1st order)
    mu = 70e6
    E1 = 70e6
    tau1 = 1e9 / E1
    E2 = 20e6
    tau2 = 1e9 / E2
    E3 = 10e6
    tau3 = 1e9 / E3
    k = 44e6

    cylinder.addObject(
        "TetrahedronViscoHyperelasticityFEMForceField",
        template="Vec3d",
        name="FEM",
        src="@topo",
        materialName="SLSNeoHookeanFirstOrder",
        ParameterSet=f"{mu} {E1} {tau1/10.0} {k}",
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
    modelVisu3.addObject("OglModel", src="@loader", color=[1.0, 0.0, 0.0, 1.0])
    modelVisu3.addObject("BarycentricMapping")

    return rootNode
