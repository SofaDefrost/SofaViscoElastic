import os
import math
import numpy as np

import Sofa.Core

# Base path dello script
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Cartella per i risultati
PLOT_DIR = os.path.join(BASE_DIR, 'plot')
os.makedirs(PLOT_DIR, exist_ok=True)


class CylinderController(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.time = 0.0
        self.node = kwargs['node']
        self.pos3 = kwargs['pos3']

        self.max1 = 0.0
        self.posmax1 = 0

        # Trova il nodo con z massima
        for j in range(len(self.pos3)):
            if self.pos3[j][2] >= self.max1:
                self.max1 = self.pos3[j][2]
                self.posmax1 = j

        print("Indice nodo max z:", self.posmax1)

        self.lin = self.node.cylinder.tetras.position.value[self.posmax1][2]

        if len(self.node.cylinder.tetras.position.value) > 4:
            epsilon = (self.node.cylinder.tetras.position.value[4][2] - self.lin) / self.lin
        else:
            epsilon = 0.0
            print("[CylinderController] ATTENZIONE: meno di 5 nodi, epsilon iniziale posto a 0.")

        with open(os.path.join(PLOT_DIR, "Burgers_new.txt"), "w") as f:
            f.write(f"{self.time} {epsilon*100}\n")

    def onAnimateBeginEvent(self, event):
        self.time = self.node.time.value

        if len(self.node.cylinder.tetras.position.value) > 4:
            epsilon = (self.node.cylinder.tetras.position.value[4][2] - self.lin) / self.lin
        else:
            epsilon = 0.0

        # Carico: step su CFF in z
        self.node.cylinder.CFF.totalForce.value = [
            0.0,
            0.0,
            (3.141592653589793e2 / 2.0) * np.heaviside(self.time, 0.0),
        ]

        print("epsilon [%] =", epsilon * 100.0)

        with open(os.path.join(PLOT_DIR, "Burgers_new.txt"), "a") as f:
            f.write(f"{self.time} {epsilon*100}\n")


def createScene(rootNode):
    rootNode.addObject(
        'VisualStyle',
        displayFlags='showVisualModels hideBehaviorModels hideCollisionModels '
                     'hideBoundingCollisionModels hideForceFields '
                     'hideInteractionForceFields hideWireframe'
    )

    rootNode.gravity = [0, 0, -9.81]
    rootNode.dt = 1e6 / 70e6
    rootNode.addObject('DefaultAnimationLoop', computeBoundingBox="0")
    rootNode.addObject('GenericConstraintSolver', tolerance=1e-24, maxIterations=1000)
    rootNode.addObject('OglSceneFrame', style='Arrows', alignment='TopRight')

    # --- Cylinder ---
    cylinder = rootNode.addChild('cylinder')

    cylinder.addObject('EulerImplicitSolver', name="Solver",
                       rayleighMass=0.0, rayleighStiffness=0.0,
                       firstOrder=True, trapezoidalScheme=False)

    cylinder.addObject('CGLinearSolver', name="ItSolver",
                       iterations="2500", tolerance="1e-30", threshold='1e-12')

    vtk_path = os.path.join(BASE_DIR, 'mesh', 'cylinder5296.vtk')
    cylinder.addObject('MeshVTKLoader', name='loader',
                       filename=vtk_path, translation=[0, 0.0, 0])

    cylinder.addObject('MechanicalObject', name='tetras', template='Vec3d', src='@loader')
    cylinder.addObject('TetrahedronSetTopologyContainer', name="topo", src='@loader')
    cylinder.addObject('TetrahedronSetTopologyModifier', name="Modifier")
    cylinder.addObject('TetrahedronSetGeometryAlgorithms', template="Vec3d", name="GeomAlgo")

    cylinder.addObject('UniformMass', totalMass="1e-6", src='@topo')

    cylinder.addObject('BoxROI', name='boxROI1',
                       box="-0.011 -0.011 -0.001  0.011 0.011 0.001", drawBoxes=True)
    cylinder.addObject('FixedProjectiveConstraint', name='fix1',
                       indices='@boxROI1.indices')

    cylinder.addObject('BoxROI', name="boxToPull",
                       box=[-0.011, -0.011, 0.05, 0.011, 0.011, 0.051],
                       drawBoxes=False)

    G1 = 70e6
    G2 = 20e6
    G3 = 10e6
    tau1 = 1e6 / G1
    tau2 = 1e6 / G2
    tau3 = 1e6 / G3
    nu = 0.44
    lamb = 1000 * G1
    lamb1 = (2* G1* (1 + nu))/(3 - 6* nu) 

    cylinder.addObject(
        'TetrahedronViscoelasticityFEMForceField',
        template='Vec3d',
        name='FEM',
        src='@topo',
        materialName="MaxwellFirstOrder",
        ParameterSet=f"{G1} {tau1} {lamb1}"
    )

    cylinder.addObject(
        'ConstantForceField',
        name="CFF",
        listening=True,
        totalForce=[0, 0, 0],
        template="Vec3d",
        src="@topo",
        indices='@boxToPull.indices'
    )

    cylinder.addObject('GenericConstraintSolver')

    # Controller
    cylinder.addObject(
        CylinderController(
            node=rootNode,
            pos3=rootNode.cylinder.tetras.position.value
        )
    )

    # Visualizzazione
    modelVisu3 = cylinder.addChild('visu')
    stl_path = os.path.join(BASE_DIR, 'mesh', 'cylinder5296.stl')
    modelVisu3.addObject('MeshSTLLoader', name='loader',
                         filename=stl_path, translation=[0.0, 0.0, 0.0])
    modelVisu3.addObject('OglModel', src='@loader', color=[1, 1, 1, 1])
    modelVisu3.addObject('BarycentricMapping')

    return rootNode
