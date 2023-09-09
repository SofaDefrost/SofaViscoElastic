# SofaViscoElastic

## Table of contents
* [Introduction](#introduction)
* [Installation](#installation)
* [Python Functions and Bindings](#python-functions-and-bindings)
* [Algorithm](#algorithm)
* [Examples](#examples)

## Introduction
SofaViscoElastic is a plugin for the Software Open Architecture Framework (SOFA) which implements the fundamental linear viscoelastic constitutive laws applied to tetrahedral meshes.
Viscoelasticity is a property of elastomeric materials that influences their mechanical behavior under dynamic conditions. In fact, viscoelastic constitutive equations are dependent on the stress/strain rate. At low stress/strain rates, a viscoelastic material behaves like a viscous liquid-like material, while at high stress/strain rates, the same material behaves like a Hookean solid. In fact, the simplest viscoelastic models are:

![Basic Models](./img/img1.png)

These two models represent the basic unit that constitutes the viscoelastic materials. They are composed of an elastic part described by the spring symbol and a viscous one represented by the dashpot.
Anyway, the Maxwell and the Kelvin-Voigt models describe the behavior of a few kinds of materials, like silly-putty, gels, etc. Furthermore, they are theoretical models that are not stable under creep (Maxwell) or stress relaxations (Kelvin-Voigt) conditions. 
Elastomers and rubbers are polymeric materials present in nature but are also used in several industrial applications. Many research fields are involved in developing and using new elastomeric materials and rubbers, such as soft robotics and surgical applications. For this reason, this plugin is indicated for users who want a realistic mechanical simulation of these materials afflicted by the viscoelastic effect.
To describe their viscoelastic properties, different viscoelastic models have to be used, like the Standard Linear Solid (SLS) Maxwell/Kelvin-Voigt representation:

![SLS Models](./img/img2.png)

They add another spring in parallel (Maxwell representation) or in series (Kelvin-Voigt representation) to make the model stable under creep and stress relaxation and are excellent for describing viscoelastic polymer rheology. In the SofaViscoElastic plugin, 9 different viscoelastic models are presented. For more theoretical information the users can refer to the paper "Considering the viscoelastic effects in soft robotic modeling" by Ferrentino et al. submitted in Soft Robotic Journal (SORO). 

## Installation
This plugin is available only for Ubuntu/Linux OS. The only dependency is the SOFA plugin "SofaPython3" (make sure it is installed).
To install this plugin from the source the user has to download this folder and place it in:
```
 $ /home/adminName/sofa/src/applications/plugins
```
Then the user has to write this in the CMakeLists.txt present in the previous destination:
```
$ sofa_add_subdirectory(plugin SofaViscoElastic SofaViscoElastic)
```
then recompile SOFA and it should start its installation. Enjoy!
For any problem contact the author at this email: pasquale.ferrentino@vub.be.

## Python Functions and Bindings
The principal function of this plugin is the so-called TethrahedronViscoelasticityFEMForceField which applies the viscoelastic constitutive law to the tetrahedral mesh uploaded in SOFA, the syntax in Python is the following :

![Python function](./img/img3.png)

The additional fields to fill in are:
* template: related to the DOF expressed in the Mechanical Object in SOFA.
* name: The name chosen by the user for the function (will appear in the SOFA simulation Graph)
* materialName: The name of the viscoelastic model that the user wants to use, he can choose between:
  - MaxwellFirstOrder
  - KelvinVoigtFirstOrder
  - SLSMaxwellFirstOrder
  - SLSKelvinVoigtFirstOrder
  - Burgers
  - MaxwellSecondOrder
  - KelvinVoigtSecondOrder
  - SLSMaxwellSecondOrder
  - SLSKelvinVoigtSecondOrder 
* ParameterSet: the lists of the material constants proper of the viscoelastic model chosen by the user. In particular, the user has to define the Young Moduli ($E_i$) and the relaxation times ($&tau;_i$) defined as the ratio between the viscosity ($&eta;_i$) and the relative Young modulus. At the end the user has to specify the Poisson Ratio (&nu;).
Furthermore in the plugin are integrated some Python Bindings to export some internal parameters of specific tetrahedrons:

![Python Binding](./img/img4.png)

The user can choose these quantities:
* getShapeVector(): Get the shape Vector of the specific tetrahedron (specified in the brackets).
* getFiberDirection(): get the fiber direction of the tetrahedron.
* getVolume(): get the volume of the tetrahedron.
* getRestVolume(): get the rest volume of the tetrahedron.
* getVolScale(): get the volume scale of the tetrahedron.
* getF(): get the deformation gradient (F) of the tetrahedron.
* getSPKStress(): get the stresses of the Second-Piola Kirchhoff tensor.
* getCauchyStress(): get the stresses Cauchy tensor.
* getVonMisesStress(): get the Von Mises stresses.

P.S. The stresses are per Element not per Node.
## Algorithm
In the Figure.4 of the paper, it is represents the SLS-Maxwell model of first order, used as example to make understand the algorithm 
## Examples
 
