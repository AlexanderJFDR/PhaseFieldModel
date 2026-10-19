# Phase Field Fracture Method (PFM) — Learning Repository

This repository documents the author's learning journey of the **Phase Field Method (PFM)** for fracture modeling. It contains a comprehensive Chinese-language summary document, numerical examples implemented in **Abaqus UEL** and **FreeFEM++**, and serves as a personal knowledge base covering the theoretical foundations, numerical implementations, and advanced topics of phase field fracture modeling.

> **⚠️ Disclaimer:** This project was compiled over a long period (starting from 2024) and may contain inconsistencies in variable notations, non-sequential reference ordering, and potential errors or omissions. Corrections and suggestions are warmly welcome.

---

## Summary Document at a Glance

| Item | Details |
| --- | --- |
| **Title** | 断裂相场法总结 (*Summary of Phase Field Fracture Method*) |
| **Author** | Heliang You, Harbin Institute of Technology |
| **Date** | October 18, 2026 |
| **Length** | 46 pages, 14 chapters, 20 figures, 47 references |
| **Language** | Chinese |
| **File** | `断裂相场法总结.pdf` |

The PDF is the current and authoritative version of the document. It extends the earlier outline to **14 chapters**, including new sections on the second-order stability of phase field solutions, advanced solution algorithms, dynamic fracture pitfalls, and length-scale-insensitive phase field models.

---

## Repository Structure

```
断裂相场法/
│
├── 断裂相场法总结.pdf            # Current summary of PFM (Chinese, 46 pages)
│
├── abaqus代码/
│   └── AT2相场拉伸/              # Abaqus UEL implementation (AT2 model)
│       ├── AT2.for               # Fortran UEL subroutine (Q4, 12 DOFs)
│       ├── origintension.inp     # Abaqus input file
│       └── origintension.odb     # Output database (example result, Git LFS)
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
            └── phasefield/       # damage*.vtu files for visualization
```

> **Note:** The `.odb` file is stored with **Git LFS** (see `.gitattributes`). Install and initialize Git LFS before cloning if you need the example database.

---

## Contents of the Summary Document (`断裂相场法总结.pdf`)

The core of this repository is the Chinese-language document **"断裂相场法总结"** (Summary of Phase Field Fracture Method), which systematically covers the following topics across **14 chapters**:

### Chapter 1: Classical Phase Field Fracture Theory (一、经典断裂相场理论)
- **1.1 AT2 Phase Field Model** — Total potential energy composed of strain energy, crack surface energy, and external work; degradation function *g(d)* and crack density function *γ(d, ∇d)*; derivation of the equilibrium and phase-field evolution equations from the first-order variational inequality with an irreversibility constraint. Both the strong form and the natural boundary conditions are given.
- **1.2 General Form of the Phase Field Equation** — Wu Jianying's unified phase-field framework with a general crack geometric function *α(d)*, scaling coefficient *c₀*, and length scale *l₀*; the displacement equation is kept in its standard variational form.

### Chapter 2: Crack Geometric Function (二、裂纹几何函数)
- **2.1 Scaling Coefficient *c₀*** — Derivation from crack-surface regularization; the 1D projection of the crack density function determines *c₀*, which affects the crack profile but not the softening law or fracture energy.
- **2.2 Crack Geometric Functions *α(d)*** — AT2 (*α = d²*), AT1 (*α = d*), and the PF-CZM interpolation (*α = ξd + (1 − ξ)d²*); 1D crack profiles and their half-bandwidths *D* (AT2: ∞, AT1: 2*l₀*, PF-CZM: *l₀*); comparison of geometric functions proposed by different researchers.
- **2.3 Singular Solution Analysis** — The homogeneous solution obtained when the strain energy vanishes; this solution provides the initial phase-field value. AT2 has the singular solution *d = 0*; AT1 has no singular solution and may exhibit numerical blow-up; PF-CZM has *d = 1*, so a thresholded history field is required to start from an undamaged state.

### Chapter 3: Degradation Function (三、退化函数)
- **3.1 Crack Initiation Strength Control** — 1D analytical derivation of the initiation strength. AT2 has zero initiation strength (damage begins immediately), whereas AT1 and PF-CZM have a finite initiation strength.
- **3.2 Ultimate Strength Control** — Stress–strain peak analysis; the damage value `d*` at peak stress; comparison of quadratic, cubic, and quartic degradation functions and the predicted ultimate strengths of AT2, AT1, and PF-CZM.
- **3.3 Damage Softening Behavior Control** — Tangent stiffness derivation and the snap-back condition; Borden's cubic degradation function *g(d; m)* and Sargado's exponential form *g(d; k, n, w)* for tuning the linear range, peak, and post-peak slope.
- **3.4 Cohesive Softening Behavior Control** — Wu's unified PF-CZM with a three-parameter degradation function; a Volterra integral-equation approach that yields an exact degradation function for a prescribed cohesive softening law; equivalent softening curves and their relation to *Gc*.
- **3.5 Phase Field Width Control** — The intrinsic material length *lch = EGc/σ²*; decoupling of the phase-field width from the ultimate strength in PF-CZM via *a₁*, limited by the positive-definiteness requirement of the iteration matrix; Lo's exponential degradation function with parameter *s*, which permits a larger phase-field width for macroscopic structures.
- **3.6 Semi-bandwidth Length** — Relation between the damage-band half-width *D* and the maximum damage value `d*`; bandwidth evolution during damage growth (increasing for PF-CZM, decreasing for Lo's model, the latter being inconsistent with damage irreversibility).
- **3.7 Stability of Solutions** — Second-order stability analysis in 1D for PF-CZM. A permissible perturbation is constructed to show that the homogeneous solution becomes unstable for sufficiently long bars; Legendre–Hadamard conditions and coercivity provide sufficient conditions for the stability of localized solutions.

### Chapter 4: Energy Decomposition Methods (四、能量分解方式)
- **4.1 Spectral Decomposition (Miehe)** — Principal-strain split into tensile and compressive parts; detailed tangent stiffness derivation including fourth-order tensor calculus; Dijk's extension to orthotropic materials.
- **4.2 Volumetric–Deviatoric Decomposition (Amor)** — Split based on volumetric (*tr ε*) and deviatoric strain; stress update and tangent stiffness; Nguyen's application to anisotropic materials.
- **4.3 Projection-based Decomposition** — Positive/negative projection of stress or strain (orthogonal decomposition); energy-norm-based projection by Wu & Cervera; orthogonality with respect to the compliance matrix; damage surfaces in principal-stress space.
- **4.4 Energy Decomposition for Anisotropic Materials** — Equivalent strain tensor via the square root of the stiffness matrix; orthogonal projection decomposition (Nguyen et al.); volumetric–deviatoric split of the equivalent strain (Ziaei-Rad et al.); properties of the projection operators and stiffness recovery.
- **4.5 Maximum History Variable** — Miehe's history field for irreversibility and phase-field bounds; modification for PF-CZM so that *d* starts from 0; Wu's history variable with a regularized crack driving force; conditions the history variable must satisfy.

### Chapter 5: Solution Methods (五、求解方法)
- **5.1 Monolithic Solution** — Fully coupled Newton–Raphson with off-diagonal coupling terms; non-convexity of the coupled system leads to convergence difficulties.
- **5.2 Staggered (Alternate Minimization) Solution** — Bourdin's approach decoupling displacement and phase field; convex sub-problems with robust convergence; multiple iterations per load increment.
- **5.3 BFGS Algorithm** — Rank-2 quasi-Newton update of the inverse stiffness matrix; initial tangent with coupling terms omitted to preserve positive definiteness; fewer stiffness factorizations at the cost of more iterations.
- **5.4 Linearized Solution** — Freezing the nonlinear degradation/geometric terms at the previous increment; full versus partial linearization and the resulting accuracy.
- **5.5 Reduced-Space Active-Set Newton Method** — Enforcing the phase-field bounds [0, 1]; definition of the active set, Newton iteration on the reduced subspace, projection onto the bounds, and active-set updates.
- **5.6 Path-Following Strategies** — Crack-surface and displacement control via a Lagrange multiplier; scaling the load so that the crack advances by a controlled amount per increment; condensed solution procedure.
- **5.7 Restrictions on the Degradation Function** — Symmetric positive definiteness condition *g″(d)H + (Gc/c₀l₀)α″(d) ≥ 0*; admissible PF-CZM length scale; Taylor-expansion modification of *g″(d)* to recover a convergent iteration matrix.
- **5.8 Iteration Termination Criteria** — Residual-based (Seleš et al.), phase-field increment-based (Bourdin), and energy-based (Ambati) convergence checks.

### Chapter 6: Adaptive Mesh Refinement (六、自适应网格划分)
- Most adaptive PFM implementations refine meshes by inserting hanging nodes, typically within AMR frameworks such as **deal.II**.
- **FreeFEM++** does not support hanging nodes; it uses the built-in `adaptmesh` function instead, which regenerates the mesh and interpolates the solution. The chapter discusses the corresponding strategy and the interpolation errors that must be controlled.

### Chapter 7: Computational Software (七、计算求解软件)
- **7.1 FreeFEM++** — Open-source FEM workflow: parameter/function definition, mesh generation, weak-form definition, staggered iteration, and post-processing; PETSc-based parallel computing; a unified framework in which AT2, AT1, and PF-CZM can be selected.
- **7.2 Abaqus Implementation** — Detailed UEL and UMAT implementation techniques.
  - **7.2.1 UMAT Subroutine** — Thermo-mechanical analogy: the phase field is stored in the temperature degree of freedom; the phase-field equation is rewritten as a nonlinear heat-conduction equation with a nonlinear heat source.
  - **7.2.2 UEL Subroutine** — Q4 and Q8 element formulations; shape functions, Jacobian, and second derivatives for Q8; displacement and phase-field element stiffness/residual; 2×2 Gaussian quadrature; degree-of-freedom reordering for the Abaqus interface.
  - **7.2.3 Practical Considerations** — Mesh size at the crack path (< *l₀*/3); linear elements (CPS3/CPS4) for phase-field convergence; variable settings (NT11 stores the phase field); convergence criteria for cohesive (PF-CZM) versus brittle (AT2) models.

### Chapter 8: Anisotropic Material Applications (八、各向异性材料的应用)
- Covers stiffness anisotropy, fracture of heterogeneous constituents, and fracture-toughness anisotropy — the most widely studied case.
- **8.1 Fourth-Order Tensor Damage Variable (Petrini)** — Fourth-order anisotropic degradation tensor *G*; direction-dependent stiffness degradation instead of isotropic damage.
- **8.2 Matrix vs. Fiber Fracture (Song)** — Dual phase-field variables for fiber and matrix failure in composite laminae; strain-energy decomposition in the principal material coordinates.
- **8.3 Modified Anisotropic Phase Field Model** — Directional crack density function containing the fiber orientation *N* and the term *β cos⟨N, ∇d⟩*; strain-energy split into fiber failure, matrix Mode I, and matrix Mode II parts; separate critical energy release rates for each mode.
- **8.4 Structural Tensor for Fracture-Toughness Anisotropy** — Second-order structural tensor *L* in the phase-field equation; derivation of the equivalent critical energy release rate *Gc(θ)* for a crack propagating at an angle *θ*; maximum/minimum fracture energies along the principal directions; extension to fourth-order tensors for strong anisotropy and zigzag crack paths.

### Chapter 9: Mixed-mode Fracture (九、断裂相场用于混合模式断裂)
- Mixed-mode fracture generally considers Mode I/II fracture (tension/shear separation), general failure criteria (tension/compression/shear), or composite failure (fiber tension, matrix shear, interface). Such problems are usually addressed by modifying the crack driving force or by using dual phase fields. Mixed-mode PFM models generally lack variational consistency, so the hybrid phase-field format is commonly used.
- **9.1 Double Phase Field Model for Rocks (Fan & Choo)** — Distinguishing tensile (Mode I) and shear (Mode II) cracks; an F-criterion determines the crack orientation; PF-CZM extended with two phase fields.
- **9.2 Mixed Power-law Phase Field Model (Shen)** — Power-law mixed-mode fracture criterion; separate history variables for tensile and shear strain energy; transition of the crack topology from Mode I to Mode II.
- **9.3 General Damage-Criterion-based Mixed-mode PFM** — Incorporation of general damage criteria (e.g., quadratic failure criteria for concrete) by replacing the equivalent stress in the crack driving force; usually implemented in the hybrid phase-field format because an explicit tangent stiffness is difficult to derive.

### Chapter 10: Dynamic Phase Field Fracture (十、动态断裂相场)
- **10.1 Implicit Dynamic** — Inertia term, initial and boundary conditions; HHT-α time integration; coupled tangent stiffness and residual for Newton–Raphson; Abaqus UEL implementation.
- **10.2 Explicit Dynamic** — Viscous regularization of the phase-field equation; central-difference update of displacement; forward-Euler update of the phase field; lumped mass and capacity matrices; Abaqus VUEL implementation with parallel computing.
- **10.3 Practical Issues in Dynamic Fracture** —
  - *Stress-wave speed mismatch*: degrading stiffness while keeping density constant lowers the wave speed near the crack tip and limits the energy available for crack growth; degrading the mass as well preserves the wave speed and reproduces crack-tip instabilities.
  - *Spurious widening of the damage zone*: without mass degradation, additional elastic waves are generated in the damaged region, causing continued fracture dissipation and broadening of the phase field.
  - *Overestimated fracture energy* (Kalthoff–Winkler example): the chapter analyzes mesh-size effects, mass degradation, irreversibility constraints, Γ-convergence, and the history variable. With these effects corrected, the integrated fracture energy is 1968 J versus a geometric value of 1808 J (an overestimate of only 8.8%).

### Chapter 11: Length-Scale-Insensitive Phase Field Models (十一、长度尺度不敏感的相场模型)
- **Motivation** — In classical AT1/AT2 models the ultimate strength depends on the phase-field width, and the admissible PF-CZM width is limited by the positive definiteness of the iteration matrix (Section 5.7). For structures much larger than the intrinsic length *lch = EGc/σ²*, this would require an excessively fine mesh. The chapter starts from the 1D decomposition of total strain into elastic and crack-induced (eigen)strain and reviews two strategies: reducing stiffness (increasing compliance) versus introducing an independent eigenstrain.
- **11.1 Classical Length-Scale-Insensitive Model** — Splits the total energy into elastic, fracture, crack, and external parts, with an additional degradation function *ω(d)* and a critical fracture energy. The 1D stress–strain response becomes linear up to a crack-initiation strain and is independent of the phase-field width. However, a snap-back can occur when *l₀ > 2Gc/(Eεc²)*, in which case arc-length control is required.
- **11.2 Improved Length-Scale-Insensitive Model** — Extends the previous model to arbitrary softening laws through two degradation functions *g₁(d)*, *g₂(d)* and critical stress/strain parameters. The initiation stress and peak stress are independent of *l₀*, and the formulation avoids the snap-back problem.
- **11.3 Eigenstrain-Based Length-Scale-Insensitive Model** — Treats the crack-induced eigenstrain *εc* as an independent, non-negative variable (analogous to plastic strain) and does not degrade the stiffness. The KKT conditions provide the evolution equations for *εc* and *d*; the stress remains elastic, the strength is independent of *l₀*, and arbitrary strength surfaces can be represented through the initiation stress *σ₀*.

### Chapter 12: Phase Field Models with Interfaces (十二、考虑界面的相场断裂模型)
- **Reserved section** — The current PDF contains only the chapter heading; content is planned for a future update.

### Chapter 13: Fatigue Phase Field Fracture Model (十三、疲劳断裂相场模型)
- **Reserved section** — The current PDF contains only the chapter heading; content is planned for a future update.

### Chapter 14: Fourth-Order Phase Field Model (十四、四阶相场模型)
- Fourth-order phase-field models introduce higher-order boundary conditions and require *C¹*-continuous elements when solved directly, so they are usually implemented with *C¹* finite elements, isogeometric analysis, or mixed reduced-order methods. Their main advantage is the representation of **strong anisotropy**.
- **14.1 Governing Equations and Boundary Conditions** — Fourth-order crack surface energy containing both *∇d ⊗ ∇d* and *Δd Δd* terms, with fourth-order structural tensors *L* and *A*; detailed derivation of the natural boundary conditions via normal/tangential decomposition of *∇δd*, Stokes' theorem, and corner conditions.
- **14.2 Strong Anisotropy** — The fourth-order structural tensor makes the equivalent fracture toughness periodic with **2π** (instead of *π* for the second-order model), allowing up to four minima. The resulting *Gc* surface can be non-convex, producing forbidden crack-propagation directions and zigzag crack patterns observed in strongly anisotropic materials.

### References (参考文献)
- 47 references covering the theoretical foundations, numerical methods, anisotropic/mixed-mode/dynamic fracture, length-scale-insensitive models, and fourth-order phase-field models.

---

## Numerical Examples

### Abaqus UEL — AT2 Phase Field Tension
- **Location:** `abaqus代码/AT2相场拉伸/`
- **Description:** A Fortran UEL subroutine implementing the AT2 isotropic phase-field model for a tension test. The element is a 4-node quadrilateral (Q4) with **12 element DOFs** (displacement *u*, *v* and phase field *d* at each node), 2×2 Gaussian integration, plane-stress stiffness input (`A11`, `A12`, `A22`, `A66`), fracture energy `Gf`, and length scale `lb`. A hybrid/staggered (alternating minimization) solution scheme is used.
- **Files:** UEL source code (`AT2.for`), Abaqus input file (`origintension.inp`), and sample results (`origintension.odb`).
- **Note:** In the current snapshot, `origintension.inp` contains the nodal data only; no element, section, boundary-condition, or step blocks were found. Complete those blocks before running the job.

### FreeFEM++ — Unified Phase Field Model Framework
- **Location:** `freefem代码/UnifiedPFM/`
- **Description:** A FreeFEM++ implementation built on the unified PFM framework. Model parameters are selected in `main.edp` through `ksi` and `c0`:
  - AT2: `ksi = 0`, `c0 = 2`
  - AT1: `ksi = 1`, `c0 = 8/3`
  - PF-CZM: `ksi = 2`, `c0 = π` (default in the script)
- **Features:**
  - Unified crack geometric function `α(d) = ξd + (1 − ξ)d²` and three-parameter degradation function
  - Energy-decomposition macros: spectral (Miehe), volumetric–deviatoric (Amor), and no decomposition (the script currently uses the VD form)
  - Staggered solution with a maximum history variable and bound enforcement of the phase field
  - Mesh templates: center-crack tension, shear, original (straight) crack, hole, and rectangular plate
  - Force–displacement output and VTK (`damage*.vtu`) output for visualization
  - PETSc sparse solvers; the script loads the `iovtk`, `gmsh`, and `ff-Ipopt` plugins
- **Key Scripts:**
  - `main.edp` — Material parameters, model selection, weak forms, staggered iteration, and output
  - `meshmodel.edp` — Mesh-generation functions for the specimen geometries

---

## Getting Started

### Prerequisites
- **Abaqus** (for the UEL example)
- **FreeFEM++** with the `iovtk`, `gmsh`, and `ff-Ipopt` plugins (for the FreeFEM example)
- **Gmsh** (only needed to regenerate `tension.msh`; the mesh is already included)
- **Git LFS** (optional, required to retrieve the example `.odb`)

### Running FreeFEM Examples
```bash
cd freefem代码/UnifiedPFM
FreeFem++ main.edp
```
Results are written to `output/displacement-force.txt` and `output/phasefield/damage*.vtu`.

### Running Abaqus UEL Examples
1. Place `AT2.for` in the working directory together with `origintension.inp`.
2. Run with Abaqus:
```bash
abaqus job=origintension user=AT2.for
```

> The UEL assembles three DOFs per node in the local order (*u*, *v*, *d*) and reorders the element residual/stiffness through the transformation matrix `QQ` before returning them to Abaqus. See Section 7.2.3 of the summary for mesh-size, element-type, and convergence recommendations.

---

## Related Repositories

- [Thermo-Mechanic Phase Field Model](https://github.com/AlexanderJFDR/Thermo-Mechanic-PhaseField-Model.git)
- [Spall Phase Field VUEL](https://github.com/AlexanderJFDR/spall_phase_field_VUEL.git)
- [PF-CZM Abaqus (Wu Jianying)](https://github.com/jianyingwu/pfczm-abaqus)

---

## Language Note

- The current summary document (`断裂相场法总结.pdf`) and all equation derivations are available **only in Chinese**.
- This README provides an English overview of the document's structure; it is not a full translation.
- The older `.docx` filename is obsolete and has been replaced by the PDF in this repository.
- An English version of the summary may be released in the future — depending on the pace of AI-assisted translation, as the author may not have time for manual translation.

---

## License

This project is intended for educational and research purposes.

---

## Acknowledgments

- Special thanks to **Prof. Wu Jianying** for providing the PF-CZM Abaqus code, which served as a reference for the UEL implementation.
- Thanks to all the researchers whose works are cited in the summary document.
