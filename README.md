# ***Crystal*** — usage guide

***Crystal*** is a deterministic lattice-traversal tool for the
global exploration of multidimensional surfaces.
Starting from a single seed point (the *root*), it systematically
walks a discrete integer lattice (grid) on which a user-defined
*cost* function is evaluated.
The traversal expands outward layer by layer and retains every
point whose cost stays below a prescribed *cost cutoff*.
Because it is exhaustive within the cutoff (rather than
stochastic), it covers all cost-relevant regions of the surface
reachable from the root.
A full description of the algorithm is given in Part I,
Sec. 2.2 (see the Citation section for Part I and Part II).

The code is organized so that the generic lattice-traversal core
is fully separated from the *system-dependent* part.
A "system" (`sys`) supplies the concrete meaning of a lattice
point: how a geometry is built from a lattice vector, how the
cost is evaluated, and what to do with the accepted points.
Systems live under `src/sys/`, and the `sys/` layer is
**arbitrarily extensible** — new systems can be plugged in
without touching the traversal kernel.
At present, a single system is shipped:

- **`rigidMC`** — general system for the treatment of **rigid
  molecular clusters**, where the cost function is the
  interaction energy (relative to the energy of the fully
  dissociated and relaxed cluster).

This document covers three topics:

1. [Compilation](#1-compilation)
2. [The `rigidMC` system](#2-the-rigidmc-system) — input files,
   instructions for running, description of hard-coded sentinel
   parameters
3. [Adding new systems](#3-adding-new-systems)

## Contact

Roland Tóbiás (rtobias@uvm.edu),
Tarun Roy (Tarun.Roy@uvm.edu), and
Bill Poirier (Bill.Poirier@uvm.edu).

## Citation

If you use ***Crystal***, please cite the following works:

- A. Pandey, B. Poirier, R. Tóbiás, T. Roy, P. Tschaikowsky,
  *****Crystal***: search engine for global exploration of
  multidimensional surfaces*, GitHub repository (2026).
- **Part I** — R. Tóbiás, T. Roy, A. Pandey, P. Tschaikowsky,
  R. Liang, B. Poirier, *Global exploration of stationary-point
  configurations in molecular clusters using the ***Crystal***
  algorithm. I. Theory and implementation*, J. Comput. Chem.
  (2026), *under review*.
- **Part II** — T. Roy, R. Tóbiás, A. Pandey, P. Tschaikowsky,
  R. Liang, B. Poirier, *Global exploration of stationary-point
  configurations in molecular clusters using the ***Crystal***
  algorithm. II. A case study of N2 clusters*, J. Comput. Chem.
  (2026), *under review*.

The ***Crystal*** algorithm itself was originally introduced in:

- A. Pandey, B. Poirier, *An algorithm to find (and plug)
  "holes" in multi-dimensional surfaces*, J. Chem. Phys.
  **152**, 214102 (2020).
  [doi:10.1063/5.0005681](https://doi.org/10.1063/5.0005681)
- A. Pandey, B. Poirier, *Plumbing potentials for molecules with
  up to tens of atoms: how to find saddle points and fix leaky
  holes*, J. Phys. Chem. Lett. **11**, 6468–6474 (2020).
  [doi:10.1021/acs.jpclett.0c01435](https://doi.org/10.1021/acs.jpclett.0c01435)
- M. Aarabi, A. Pandey, B. Poirier, *"On-the-fly" Crystal: How
  to reliably and automatically characterize and construct
  potential energy surfaces*, J. Comput. Chem. **45**,
  1261–1278 (2024).
  [doi:10.1002/jcc.27324](https://doi.org/10.1002/jcc.27324)

---

## 1. Compilation

Compilation is driven by the `Makefile` in the root directory.

```sh
make          # build crystal.x
make clean    # remove the contents of build/ and crystal.x
```

Behavior is controlled by two variables at the top of the
`Makefile`:

| Variable | Allowed values | Meaning |
|----------|----------------|---------|
| `FC`     | `ifx` (default), `ifort`, `gfortran` | Fortran compiler |
| `MODE`   | `release` (default), `debug` | optimized vs. checked/debug build |

These can be overridden from the command line, e.g.:

```sh
make FC=gfortran MODE=debug
```

### Prerequisites

- **Fortran compiler:** Intel oneAPI (`ifx`/`ifort`) or GNU
  (`gfortran`).
  Tested with Intel oneAPI 2024.2.1 (`ifort` and `ifx`) and
  2025.3 (`ifx`), and with `gfortran` 11.5.0, 12.4.0 and 14.2.0.
  `gfortran` 10 and earlier are not supported.
- **LAPACK/BLAS:** static libraries under
  `~/.local/lib/lapack-3.12.0/lib/$(TC)`, where `$(TC)` is the
  toolchain (`intel` or `gnu`).
  If located elsewhere, edit the `LAPACK_PATH` variable in the
  `Makefile`.
Object files (`.o`, `.mod`) go to the `build/` directory; the
final binary is `crystal.x` in the root directory.

> Inter-module dependencies are described by the `.mk` files
> under `mk/` (`obgen.mk`, `util.mk`, `depsys.mk`,
> `generic.mk`).
> When adding a new source file, its dependencies must be
> activated there — see [Section 3](#3-adding-new-systems).

<!-- -->

> **⚠️ Warning — reproducibility across builds**
>
> The lattice traversal is not guaranteed to be reproducible from
> one build to another.
> Two of its decisions are exact comparisons on computed costs:
>
> - a lattice point is retained only while its cost stays below
>   `cut` (`cost < cut`);
> - a lattice point is assigned to a cell by rounding
>   `cost/contrac` to the nearest integer.
>
> Floating-point operations are processor- and compiler-dependent,
> so a cost falling essentially exactly on either threshold may
> land on either side.
> Runs made with different builds may therefore traverse slightly
> different portions of the lattice.

---

## 2. The `rigidMC` system

`rigidMC` is the system for locating stationary points on the
interaction potential energy surfaces of rigid molecular
clusters.
Each monomer is described by 6 curvilinear rigid-body
coordinates (3 rotational Euler angles + 3 spherical
translational coordinates), and the traversal explores the
intermonomer configuration space.
The cost function is the interaction energy (relative to the
energy of the fully dissociated and relaxed cluster).

### Running

```sh
./crystal.x inp.dat > out.dat
```

The single command-line argument is the input file name (plus
the associated file path).
The number of threads is set by the `nthreads` key in the input
(OpenMP).

Ready-made examples live under `test/clus/*/` (`inp.dat`,
reference `out.dat`, `stac.xyz`), and `test/slurm/*.sh` provides
SLURM job templates.

### Structure of the input file

The input is a simple key–value format organized into sections
(blocks).
A typical `test/clus/Ne3`-style input:

```
Generic parameters:

  cut = 0.0
  contrac = 0.0
  sys_name = rigidMC
  nthreads = 60
  plot = 0

System-dependent parameters:

    Geometry:

      cluster  = /path/to/xyz/Nen/Ne3_GM.xyz
      monomer1 = /path/to/xyz/Nen/Ne.xyz
      monomer2 = /path/to/xyz/Nen/Ne.xyz
      monomer3 = /path/to/xyz/Nen/Ne.xyz

    Customize:

      custom_drive = normal
      ang_step[deg] = 115.0
      rad_step[A] = relax

    Cost:

      cost_label = Energy
      cost_unit  = cm-1
      cost_spec  = /path/to/par/Ne2.f2b
      cost_drive = f2b_homo

    Task:

      task_drive = msopt
      stac_point_file = stac.xyz
      max_targets = none
      max_ener = none
      max_order = none
      skip_flat = 0
```

#### `Generic parameters` — generic (system-independent) settings

These are read by the traversal core and apply to every system.

| Key        | &nbsp;&nbsp;&nbsp;Type&nbsp;&nbsp;&nbsp; | Meaning |
|------------|--------|---------|
| `cut`      | real   | cost cutoff (cm⁻¹); during traversal, points above `cut` are not retained |
| `contrac`  | real   | if nonzero, points are binned into cells of this width (cm⁻¹) and points falling into the same cell are merged |
| `sys_name` | string | name of the system; selects the driver set (here `rigidMC`) |
| `nthreads` | int    | number of OpenMP threads |
| `plot`     | int    | generate projection plots (0 = off) |

#### `System-dependent parameters` — `rigidMC` blocks

The `rigidMC` system expects four blocks.

**`Geometry`** — the cluster and its monomers

| Key        | &nbsp;&nbsp;&nbsp;Type&nbsp;&nbsp;&nbsp; | Meaning |
|------------|--------|---------|
| `cluster`  | string&nbsp;/&nbsp;`none` | XYZ file of the initial (root) cluster geometry; if `none`, ***Crystal*** generates a starting geometry itself for the root (details in Part I, Sec. 2.3) |
| `monomerN` | string | XYZ file of the *N*-th monomer (`monomer1`, `monomer2`, …); the number of monomers is inferred from these keys |

When `cluster` is set to `none`, no explicit root geometry is
read from disk.
Instead, ***Crystal*** generates a starting configuration on its
own and uses it as the root.
The details of how this initial geometry is produced are
described in Part I, Sec. 2.3.

**`Customize`** — lattice construction / propagation

| Key             | &nbsp;&nbsp;&nbsp;Type&nbsp;&nbsp;&nbsp; | Meaning |
|-----------------|--------|---------|
| `custom_drive`  | string | driver for the `customize` substep (currently, `normal` only) |
| `ang_step[deg]` | real   | angular lattice spacing, in ° |
| `rad_step[A]`   | real&nbsp;/&nbsp;string | radial lattice spacing, in Å; the keyword `relax` requests local radial relaxation instead of a fixed step |

**`Cost`** — the cost function (interaction energy for `rigidMC`)

| Key          | &nbsp;&nbsp;&nbsp;Type&nbsp;&nbsp;&nbsp; | Meaning |
|--------------|--------|---------|
| `cost_label` | string | free-text label (e.g. `Energy`) |
| `cost_unit`  | string | unit of the cost (only `cm-1` is currently supported) |
| `cost_spec`  | string | path to the potential parameter file (`.f2b`, see below) |
| `cost_drive` | string | cost-evaluation driver (only `f2b_homo` is currently implemented) |

**`Task`** — post-processing of the accepted points

| Key               | &nbsp;&nbsp;&nbsp;Type&nbsp;&nbsp;&nbsp; | Meaning |
|-------------------|--------|---------|
| `task_drive`      | string | task driver (`msopt` = multi-start optimization) |
| `stac_point_file` | string | output XYZ file of the located stationary points |
| `max_targets`     | int&nbsp;/&nbsp;`none` | upper limit on the number of points to be optimized (`none` = no limit) |
| `max_ener`        | real&nbsp;/&nbsp;`none` | upper limit on the energy of the points to be optimized (`none` = no limit) |
| `max_order`       | int&nbsp;/&nbsp;`none` | upper limit on the saddle order of the points optimized (`none` = no limit) |
| `skip_flat`       | int    | skip flat points (stationary geometries with numerically unstable Hessian eigenvalues); allowed values: 0/1 |

> The bracketed suffix `[deg]` / `[A]` is the unit of the key;
> no other units are currently supported.

### Potential parameter file (`.f2b`)

The `.f2b` files hold the parameters of the two-body
potential used by `cost_drive = f2b_homo`.
The F2B ("flexible two-body") format is a restructured
representation of the PES parametrization used in **autoPES**
(see <https://www.physics.udel.edu/~szalewic/SAPT/download.html>).

`cost_spec` points to the `.f2b` file, and `cost_drive` selects
how it is evaluated.
Currently only `f2b_homo` is implemented: a homogeneous cluster
in which the same potential-energy parametrization is used for
all monomer pairs.

### The sentinel parameters

`src/sys/rigidMC/engine/sentinel/type_sentinel_rigidMC_inc.f90`
holds the numerical tolerances and limits of the local
optimization / geometry handling **hard-coded into the source**.
These do not come from the input — changing them requires
editing the file and recompiling.
The tables below summarize what each one means (for the
underlying algorithm, see Sec. 2.3 in Part I).
These parameters govern the filtering of redundant stationary
points, the handling of the cluster geometries, and the
convergence of the local optimizations.

#### Integer parameters

| Parameter     | &nbsp;Value&nbsp; | Meaning |
|---------------|:-----:|---------|
| `mreconds`    | 200   | maximum reconditioning cycles per optimization |
| `mdcbis`      | 40    | maximum bisection steps when a trial displacement has to be scaled back |
| `msharps`     | 5     | maximum consecutive cycles without gradient improvement |
| `mbfgs_relax` | 10    | maximum BFGS iterations per monomer pair in the radial relaxation |
| `mrest_relax` | 3     | maximum restarts of the radial relaxation |
| `mbfgs`       | 10    | maximum BFGS iterations per reconditioning cycle |
| `mnp`         | 10    | maximum Newton-Powell iterations per reconditioning cycle |

#### Real parameters

| Parameter        | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Value&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Unit&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Meaning |
|------------------|:------------------:|:-------:|---------|
| `ang_f2doff`     |      `1.e-3`       | °       | angular offset for the finite-difference derivatives |
| `rad_f2doff`     |      `1.e-5`       | Å       | radial offset for the finite-difference derivatives |
| `rad_gtol`       |       `0.05`       | cm⁻¹/Å  | radial gradient tolerance |
| `ang_gtol`       |      `1.e-3`       | cm⁻¹/°  | angular gradient tolerance |
| `rad_dvtol`      |      `1.e-3`       | Å       | radial displacement tolerance |
| `ang_dvtol`      |      `0.01`        | °       | angular displacement tolerance |
| `ang_scal`       |       `30.0`       | °       | angular coordinate scale passed to the optimizer |
| `rad_scal`       |       `0.5`        | Å       | radial coordinate scale passed to the optimizer |
| `sing_tol`       |      `1.e-12`      | —       | tolerance for numerical singularity |
| `rig_tol`        |      `1.e-5`       | Å       | rigidity tolerance (RMSD) for checking rigid-monomer consistency |
| `flat_fac`       |       `2.0`        | —       | flatness factor |
| `collin_tol`     |      `1.e-3`       | Å²      | eigenvalue threshold for determining the shape rank (collinear / planar / 3D) of the atomic arrangement |
| `pgrp_tol`       |      `5.e-2`       | Å       | tolerance for point-group determination |
| `rsmdet_tol`     |       `0.9`        | —       | relative tolerance for selecting the anchors based on the spherical-coordinate metric tensor |
| `same_ener_rtol` |      `1.e-3`       | —       | relative tolerance for treating two stationary points as isoenergetic |
| `grmsdLB_tol`    |      `0.025`       | Å       | RMSD lower bound below which two isoenergetic stationary points are taken to be the same |
| `lam_cov`        |       `1.5`        | —       | covalent-radius scaling factor for detecting too close atoms |
| `lam_vdw`        |       `1.5`        | —       | van der Waals-radius scaling factor for detecting cluster dissociation |
| `recond_tol`     |       `0.05`       | —       | reconditioning tolerance (relative) |

> **NOTE — the `?` suffix on point-group names**
>
> A point group is accepted only once the operations found close
> into a group of the expected size.
> If they do not, the detection is repeated with a progressively
> tighter tolerance, and the name that finally closes is marked
> with a trailing `?` (e.g. `Cs?` in `stac.xyz`).
> For distorted structures, starting from a larger `pgrp_tol` may
> let the intended group be found in a single pass.

#### Derived (allocatable) fields

These are filled internally, do not edit them directly:

| Field            | Derived from | Meaning |
|------------------|--------------|---------|
| `f2doffsets(:)`  | `ang_f2doff`, `rad_f2doff` | 6 finite-difference offsets per monomer (5 angular + 1 radial) |

---

## 3. Adding new systems

The `sys/` layer is where ***Crystal*** is extended.
The generic lattice-traversal core (`src/core/`) is agnostic to
the physical problem: it communicates with a system only through
a fixed set of generic procedure pointers declared in
`src/core/mod_sys.f90`.
To add a new system `foo` alongside `rigidMC`:

### 3.1 Source files

Create a `src/sys/foo/` tree (modeled on `src/sys/rigidMC/`).
The minimum required drivers, expected by `mod_sys.f90` *via* its
abstract interfaces, are five:

| Driver procedure | Responsibility | Interface |
|------------------|----------------|-----------|
| `load_foo`       | read the system-dependent input, initialize the workspace | `pload` |
| `globalize_foo`  | map the working coordinates onto the global lattice variables | `pglobalize` |
| `customize_foo`  | build geometries / constrain the propagation | `pcustomize` |
| `cost_foo`       | evaluate the cost | `pcost` |
| `task_foo`       | arbitrary process on the accepted points (print, local optimizations) | `ptask` |

The argument lists of these subroutines must match the
corresponding abstract interfaces (`pload`, `pglobalize`,
`pcustomize`, `pcost`, `ptask`) in `mod_sys.f90`.

### 3.2 Activation in `mod_sys.f90`

A system becomes active when `init_sys()` binds the generic
procedure pointers to its drivers, keyed on `sys_name`.
Add a new branch to the `select case(sys_name)` construct:

```fortran
use mod_user_foo, only : &
   load_foo, globalize_foo, customize_foo, cost_foo, task_foo
...
select case(sys_name)
   case('rigidMC')
      ...
   case('foo')                 ! <-- new branch
      load      => load_foo
      globalize => globalize_foo
      customize => customize_foo
      cost      => cost_foo
      task      => task_foo
   case default
      call catch_error(err=.true., msg='Unknown "sys_name".', ...)
end select
```

### 3.3 Build dependencies (`mk/`)

- Create `mk/sys/foo.mk` (modeled on `mk/sys/rigidMC.mk`),
  describing *via* `obgen` the compile dependencies of the new
  modules (`mod_engine_foo.f90`, `mod_user_foo.f90`, the
  `*_inc.f90` include files).
- In `mk/depsys.mk`, add `mod_user_foo.o` to the `DEP` list of
  `mod_sys.f90`, and pull in the new rules with
  `include mk/sys/foo.mk`.

### 3.4 Running the new system

In the `Generic` block of the input, set `sys_name = foo`.
The `System-dependent parameters` section is parsed entirely by
`load_foo`, so its contents can follow any format the new system
chooses to read.
The generic keys (`cut`, `contrac`, `nthreads`, `plot`) are
unchanged.
