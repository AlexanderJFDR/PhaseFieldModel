# Phase Field Fracture Method (PFM) — Learning Repository

This repository documents the author's early-stage learning journey of the **Phase Field Method (PFM)** for fracture modeling. It contains a comprehensive summary document (in Chinese), numerical examples implemented in **Abaqus UEL** and **FreeFEM++**, and serves as a personal knowledge base for understanding the theoretical foundations, numerical implementations, and advanced topics in phase field fracture modeling.

> **⚠️ Disclaimer:** This project was compiled over a long period (starting from 2024), and may contain inconsistencies in variable notations, non-sequential reference ordering, and potential errors or omissions. Corrections and suggestions are warmly welcome.

---

## Repository Structure

```
断裂相场法/
│
├── 断裂相场法总结.docx          # Comprehensive summary of PFM (in Chinese)
│
├── abaqus代码/
│   └── AT2相场拉伸/              # Abaqus UEL implementation (AT2 model)
│       ├── AT2.for               # Fortran UEL subroutine
│       ├── origintension.inp     # Abaqus input file
│       └── origintension.odb     # Output database (example result)
│
└── freefem代码/
    └── UnifiedPFM/               # FreeFEM++ implementation (unified PFM framework)
        ├── main.edp              # Main solver script
        ├── meshmodel.edp         # Mesh generation functions
        ├── mesh/                 # Mesh files
        │   ├── tension.geo       # Gmsh geometry for tension test
        │   └── tension.msh       # Generated mesh
        └── output/               # Simulation output
            ├── displacement-force.txt
            └── phasefield/       # VTU files for phase field visualization
```

---

## Contents of the Summary Document (`断裂相场法总结.docx`)

The core of this repository is the Chinese-language document **"断裂相场法总结"** (Summary of Phase Field Fracture Method), which systematically covers the following topics across 13 major chapters:

### Chapter 1: Classical Phase Field Fracture Theory (一、经典断裂相场理论)
- **1.1 AT2 Phase Field Model** — Variational formulation based on total potential energy (strain energy + crack surface energy + external work); degradation function *g(d)* and crack density function *γ(d)*; derivation of equilibrium equation and phase field evolution equation via first-order variational inequality with irreversibility constraints.
- **1.2 General Form of Phase Field Equation** — Wu Jianying's unified phase field framework with general crack geometric function *α(d)*.

### Chapter 2: Crack Geometric Function (二、裂纹几何函数)
- **2.1 Scaling Coefficient *c₀*** — Derivation from crack surface regularization under 1D projection.
- **2.2 Crack Geometric Functions *α(d)*** — AT2 (*α = d²*), AT1 (*α = d*), and cohesive zone models; 1D crack shape profiles and semi-bandwidth lengths; comparative analysis of different geometric functions proposed by various researchers.
- **2.3 Singular Solution Analysis** — Initial homogeneous phase field solution (when strain energy = 0); stability issues for AT1 and PF-CZM.

### Chapter 3: Degradation Function (三、退化函数)
- **3.1 Crack Initiation Strength Control** — 1D analytical derivation; AT2 has zero initiation strength (damage from the start); AT1 has non-zero initiation strength.
- **3.2 Ultimate Strength Control** — Stress-strain curve peak analysis; quadratic vs. cubic vs. quartic degradation functions; comparison across AT2, AT1, PF-CZM models.
- **3.3 Damage Softening Behavior Control** — Tangent stiffness derivation; snap-back phenomenon; Borden's cubic degradation function with parameter *m*.
- **3.4 Cohesive Softening Behavior Control** — Wu's unified PF-CZM with three-parameter degradation function; Volterra integral equation-based cohesive model (Feng et al.) providing exact degradation function from prescribed softening laws.
- **3.5 Phase Field Width Control** — Length-scale insensitivity; Wu's PF-CZM achieving decoupling via parameter *a₁*; Lo's exponential degradation function enabling larger crack bandwidth for macro-structures.
- **3.6 Semi-bandwidth Length** — Relationship between damage band half-width and maximum damage value; bandwidth evolution during damage progression (increasing for PF-CZM, decreasing for Lo's model).

### Chapter 4: Energy Decomposition Methods (四、能量分解方式)
- **4.1 Spectral Decomposition (Miehe)** — Principal strain-based split into tensile/compressive parts; detailed tangent stiffness derivation including fourth-order tensor calculus; Dijk's extension to orthotropic materials.
- **4.2 Volumetric-Deviatoric Decomposition (Amor)** — Split based on volumetric (*tr(ε)*) and deviatoric strain; stress update and tangent stiffness; Nguyen's application to anisotropic materials.
- **4.3 Projection-based Decomposition** — Positive/negative stress tensor projection (orthogonal decomposition); energy-norm-based projection by Wu Jianying & Cervera; orthogonality conditions w.r.t. compliance matrix; damage surface in principal stress space.
- **4.4 Energy Decomposition for Anisotropic Materials** — Equivalent strain tensor method via stiffness matrix square root; orthogonal projection decomposition (Nguyen et al.); volumetric-deviatoric split on equivalent strain (Ziaei-Rad et al.); projection operator properties and stiffness recovery.
- **4.5 Maximum History Variable** — Miehe's history field for irreversibility and phase field bounds; modification for PF-CZM to ensure *d* starts from 0; Wu's history variable with regularized crack driving force.

### Chapter 5: Solution Methods (五、求解方法)
- **5.1 Monolithic Solution** — Fully coupled Newton-Raphson with off-diagonal coupling terms; non-convexity leading to convergence difficulties.
- **5.2 Staggered (Alternate Minimization) Solution** — Bourdin's approach decoupling displacement and phase field; convex sub-problems with robust convergence; multi-iteration per load increment.
- **5.3 BFGS Algorithm** — Rank-2 quasi-Newton update; inverse BFGS formulation; initial tangent with coupling terms omitted for positive definiteness.
- **5.4 Linearized Solution** — Freezing nonlinear terms (degradation/geometric functions); partial vs. full linearization; accuracy considerations.
- **5.5 Reduced Space Active Set Newton Method** — Handling phase field bound constraints [0,1]; active set definition and iterative projection.
- **5.6 Path-following Strategies** — Crack surface control and displacement control via Lagrange multiplier; condensed solution procedure to reduce system size.
- **5.7 Restrictions on Degradation Functions** — Condition for symmetric positive definite tangent matrix; Wu's length-scale criterion for PF-CZM; Taylor expansion modification for second derivative of degradation function.
- **5.8 Iteration Termination Criteria** — Residual-based (Karlo), increment-based (Bourdin), and energy-based (Ambati) convergence checks.

### Chapter 6: Adaptive Mesh Refinement (六、自适应网格划分)
- Hanging node-based refinement strategies; avoiding node position changes to prevent error accumulation from interpolation.

### Chapter 7: Computational Software (七、计算求解软件)
- **7.1 FreeFEM++** — Open-source FEM implementation; weak form definition; staggered iteration with PETSc parallel computing; configurable model framework (AT2/AT1/PF-CZM).
- **7.2 Abaqus Implementation** — Detailed UEL and UMAT implementation techniques.
  - **7.2.1 UMAT Subroutine** — Thermo-mechanical analogy: mapping phase field to temperature field; nonlinear heat source function.
  - **7.2.2 UEL Subroutine** — Q4 and Q8 element formulations; shape functions and Jacobian; displacement and phase field element stiffness/residual; 2×2 Gaussian quadrature; degree-of-freedom reordering for Abaqus interface.
  - **7.2.3 Practical Considerations** — Mesh size (< *l₀*/3); linear elements (CPS3/CPS4) for phase field convergence; variable settings (NT11 → phase field); convergence criteria for brittle vs. cohesive fracture.

### Chapter 8: Anisotropic Material Applications (八、各向异性材料的应用)
- Covers stiffness anisotropy, heterogeneous component fracture, and fracture toughness anisotropy — the most widely studied case.
- **8.1 Fourth-order Tensor Damage Variables (Petrini)** — Anisotropic degradation tensor *G*; direction-dependent stiffness degradation avoiding isotropic damage.
- **8.2 Matrix vs. Fiber Fracture (Song)** — Dual phase field variables for fiber and matrix failure in composite laminae; strain energy decomposition in principal material coordinates.
- **8.3 Modified Anisotropic Phase Field Model** — Directional crack density function with fiber orientation; energy decomposition for fiber failure, Mode I and Mode II matrix failure; separate critical energy release rates.
- **8.4 Directional Tensor for Fracture Toughness Anisotropy** — Second-order directional tensor in the phase field equation to represent orientation-dependent fracture toughness; derivation of equivalent critical energy release rate *Gc(θ)* as a function of crack propagation direction; extension to fourth-order tensors for strong anisotropy and zigzag crack patterns.

### Chapter 9: Mixed-mode Fracture (九、断裂相场用于混合模式断裂)
- Mixed-mode fracture generally considers Mode I/II fracture (tension/shear separation), general failure criteria (tension/compression/shear), or composite failure (fiber tension, matrix shear, interface). These problems are typically addressed by modifying crack driving forces or using dual phase fields. Notably, mixed-mode PFM models generally lack variational consistency, so the hybrid phase field approach is commonly used.
- **9.1 Double Phase Field Model for Rocks (Fan)** — Distinguishing tensile (Mode I) and shear (Mode II) cracks; F-criterion-based crack orientation determination; PF-CZM extension with dual phase fields.
- **9.2 Mixed Power-law Phase Field Model (Shen)** — Power-law mixed-mode fracture criterion; history variables for tensile and shear strain energy; crack topology transition from Mode I to Mode II.
- **9.3 General Damage Criterion-based Mixed-mode PFM** — Incorporating general damage criteria (e.g., quadratic failure criteria for concrete) into the phase field framework by replacing the equivalent stress in the crack driving force; typically implemented using the hybrid phase field format due to the difficulty of deriving an explicit tangent stiffness.

### Chapter 10: Dynamic Phase Field Fracture (十、动态断裂相场)
- **10.1 Implicit Dynamic** — Inertia term incorporation; HHT-α time integration; fully coupled tangent stiffness matrix; Abaqus UEL implementation.
- **10.2 Explicit Dynamic** — Viscous regularization of phase field; central difference for displacement; forward Euler for phase field; lumped mass/capacity matrices; Abaqus VUEL for parallel computing.

### Chapter 11: Phase Field with Interfaces (十一、考虑界面的相场断裂模型)
- Reserved section with placeholder content.

### Chapter 12: Fatigue Phase Field Fracture Model (十二、疲劳断裂相场模型)
- Reserved section with placeholder content.

### Chapter 13: Fourth-order Phase Field Model (十三、四阶相场模型)
- Reserved section with placeholder content.

### References (参考文献)

---

## Numerical Examples

### Abaqus UEL — AT2 Phase Field Tension
- **Location:** `abaqus代码/AT2相场拉伸/`
- **Description:** A Fortran UEL subroutine implementing the AT2 isotropic phase field model for a tensile test. The element formulation uses a 4-node quadrilateral (Q4) with alternating minimization, supporting hybrid phase field formulation.
- **Files:** UEL source code (`AT2.for`), Abaqus input file (`origintension.inp`), and sample results (`origintension.odb`).

### FreeFEM++ — Unified Phase Field Model Framework
- **Location:** `freefem代码/UnifiedPFM/`
- **Description:** A comprehensive FreeFEM++ implementation supporting multiple phase field models (AT2, AT1, PF-CZM) within a unified framework. Features include:
  - Energy decomposition: spectral, volumetric-deviatoric, or no decomposition
  - Staggered solution scheme with history variable
  - Multiple mesh generation templates (tension, shear, hole, plate)
  - VTK output for phase field visualization
- **Key Scripts:**
  - `main.edp` — Main solver with configurable material parameters, crack geometry functions, and degradation functions
  - `meshmodel.edp` — Mesh generation functions for various specimen geometries

---

## Getting Started

### Prerequisites
- **Abaqus** (for UEL examples)
- **FreeFEM++** (for FreeFEM examples)
- **Gmsh** (for mesh generation, optional)

### Running FreeFEM Examples
```bash
cd freefem代码/UnifiedPFM
FreeFem++ main.edp
```

### Running Abaqus UEL Examples
1. Place `AT2.for` in the working directory
2. Run with Abaqus:
```bash
abaqus job=origintension user=AT2.for
```

---

## Related Repositories

- [Thermo-Mechanic Phase Field Model](https://github.com/AlexanderJFDR/Thermo-Mechanic-PhaseField-Model.git)
- [Spall Phase Field VUEL](https://github.com/AlexanderJFDR/spall_phase_field_VUEL.git)
- [PF-CZM Abaqus (Wu Jianying)](https://github.com/jianyingwu/pfczm-abaqus)

---

## Language Note

- The summary document (`断裂相场法总结.docx`) is currently available **only in Chinese**.
- An English version may be released in the future — depending on the pace of AI-assisted translation, as the author may not have time for manual translation.

---

## License

This project is intended for educational and research purposes.

---

## Acknowledgments

- Special thanks to **Prof. Wu Jianying** for providing the PF-CZM Abaqus code, which served as a reference for the UEL implementation.
- Thanks to all the researchers whose works are cited in the summary document.
