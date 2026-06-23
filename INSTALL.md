# Legacy — Installation and Usage Guide

## What is Legacy?

**Legacy** (also known as ResBos-Legacy or Legacy Fortran) is a Fortran/C++ code that computes
resummed transverse-momentum distributions of electroweak vector bosons (W±, Z0, H0, γ, etc.)
and their associated processes in hadronic collisions. It implements the
**Collins–Soper–Sterman (CSS)** resummation formalism, combining:

- A resummed **CSS piece** (dominant at low qT)
- A fixed-order **perturbative Y-piece** (dominant at high qT)
- A **non-perturbative Sudakov form factor** (g1, g2, g3 parameters)

The code computes `d²σ / dqT dy` (and optionally integrated cross sections) on a grid of
transverse momentum (qT), rapidity (y), and invariant mass (Q) points.

**Version:** `main.for` v4.4.3 / `res.for` v5.2 / `pert.for` v6.1

---

## Code Structure

```
legacy_final_vesion/
├── main.for          # Main driver program
├── res.for           # CSS resummation kernel (v5.2)
├── pert.for          # Perturbative Y-piece (v6.1)
├── pert_vj.for       # Perturbative piece for vector+jet processes
├── pda.for           # Phase-space integration helpers
├── pion.for          # Pion PDF handling
├── common.for        # Shared COMMON blocks
├── CT14Pdf.f         # CT14 built-in PDF interface
├── lhapdf.cpp        # C++ bridge to LHAPDF6
├── EvlPac02b.for     # PDF DGLAP evolution package
├── PrzPac02b.for     # Parton distribution parametrizations
├── QcdPac02b.for     # QCD running coupling / utilities
├── UtlPac02b.for     # General utilities
├── EwkPac02b.for     # Electroweak parameters
├── Makefile          # Build system
├── legacy.in         # Main input control file  <-- EDIT THIS
├── inp/
│   ├── q_grid.inp    # Grid in invariant mass Q [GeV]
│   ├── qt_grid.inp   # Grid in transverse momentum qT [GeV]
│   └── y_grid.inp    # Grid in rapidity y
├── hoppet/           # Bundled HOPPET library source
├── CT18ResKF.00.pds  # CT18 PDF data file
├── CT18ResKF00.pds   # CT18 PDF data file (alternate)
├── pdf00.pds         # Default internal PDF data file
└── keep_codes/       # Original reference copies
```

---

## Dependencies

| Library | Version | Purpose |
|---------|---------|---------|
| gfortran | ≥ 4.8 | Fortran compiler |
| g++ | ≥ 7 (C++17) | C++ compiler for LHAPDF bridge |
| LHAPDF | 6.x | External PDF sets (CT18NNLO, etc.) |
| HOPPET | 1.x | PDF evolution (bundled in `hoppet/`) |

---

## Installation on MSU HPCC

### 1. Load required modules

Log in to HPCC and load the necessary compiler modules:

```bash
module purge
module load GCC/12.3.0       # or: module load GNU/9.3.0
module load OpenMPI/4.1.5    # if MPI is needed by other codes
```

Check available modules with `module avail gcc` or `module spider gcc`.

### 2. Install LHAPDF6

```bash
# Download source
mkdir -p /mnt/home/lopezels/SourceCodes
cd /mnt/home/lopezels/SourceCodes

wget https://lhapdf.hepforge.org/downloads/?f=LHAPDF-6.5.4.tar.gz -O LHAPDF-6.5.4.tar.gz
tar -xzf LHAPDF-6.5.4.tar.gz
cd LHAPDF-6.5.4

# Configure and install to your InstallSources directory
./configure --prefix=/mnt/home/lopezels/InstallSources/LHAPDF \
            CXX=g++ CXXFLAGS="-O2 -std=c++17"
make -j4
make install

# Download the PDF sets you need
/mnt/home/lopezels/InstallSources/LHAPDF/bin/lhapdf install CT18NNLO
/mnt/home/lopezels/InstallSources/LHAPDF/bin/lhapdf install CT14nnlo
```

### 3. Build HOPPET (bundled)

HOPPET is included in the `hoppet/` subdirectory. Build it in-place:

```bash
cd /path/to/legacy_final_vesion/hoppet

./configure --prefix=$(pwd) FC=gfortran
make -j4
make install
```

This places `libhoppet_v1.a` in `hoppet/src/`, which the Makefile expects at `-L./hoppet/src/`.

### 4. Build Legacy

Make sure your environment has LHAPDF on the path (or set the variables manually):

```bash
export PATH=/mnt/home/lopezels/InstallSources/LHAPDF/bin:$PATH
export LD_LIBRARY_PATH=/mnt/home/lopezels/InstallSources/LHAPDF/lib:$LD_LIBRARY_PATH
```

Then compile:

```bash
cd /path/to/legacy_final_vesion
make clean
make
```

A successful build produces the executable `main` in the same directory.

To verify:
```bash
./main --help   # (will just run with defaults; not a real --help flag)
ls -lh main     # should be ~3-4 MB
```

---

## Files to Copy to HPCC

You need these files and directories present together when running:

```
legacy_final_vesion/
├── main              # compiled executable
├── legacy.in         # input control file
├── inp/
│   ├── q_grid.inp
│   ├── qt_grid.inp
│   └── y_grid.inp
├── CT18ResKF.00.pds  # if using CT18 internal PDFs
├── CT18ResKF00.pds
└── pdf00.pds         # if using internal CTEQ6 PDFs
```

The LHAPDF external sets (e.g. `CT18NNLO`) are accessed from the LHAPDF system installation
directory (`/mnt/home/lopezels/InstallSources/LHAPDF/`), not from the run directory.

---

## Running Legacy Interactively

```bash
# From the run directory:
./main               # reads legacy.in, writes legacy.out

./main my_run_name   # writes output to my_run_name.out
```

The program reads `legacy.in` from the **current working directory** at startup.
Grid files referenced inside `legacy.in` (e.g. `./inp/q_grid.inp`) are relative to
the current directory as well.

---

## Understanding `legacy.in`

The first 21 lines of `legacy.in` are the active input. Lines after the `########` marker
are documentation only.

```
1,0.d0,1,0         > iBeam=-2/-1/0/1:PiN/ppB/pN/pp, FRACT_N, KinCorr, ResumType
8000.0,3,0         > ECM [GeV], LTO, iFast
172., 80.358, 91.118, 125, 0   > top/W/Z/Higgs/photon masses [GeV]
W+                 > Boson type: W+, W-, Z0, H0, A0, AA, AG, ZZ, ZG, ...
50,2,1,2           > ISTR (PDF set index), iPionPDF, i_RunMass, i_Model
CSS                > Resummation scheme: CSS or CFG
5,99,1,1,1         > iNONPERT, IFLAG_C3, C1/B0, C2, C3/B0
3,2,1,0,2          > Sudakov orders: A, B, CFns, i_FSR, DelSigOrder
1,1,1              > iscale, muR/mu, muF/mu
lha_CT18NNLO       > PDF set name (prefix lha_ means use LHAPDF6)
0.06983d0,0.06984d0,0.15d0,1.3d0,4.75d0  > u,d,s,c,b quark masses [GeV]
0.5                > bMax (b-space cutoff)
0                  > iProc=0/1/2: all/qqbar/qG
./inp/q_grid.inp   > Q (invariant mass) grid
./inp/qt_grid.inp  > qT grid
./inp/y_grid.inp   > rapidity grid
-                  > separator
0.21d0,0.68d0,-0.126d0,1.6d0,1  > g1,g2,g3,Q0,nG (BLNY nonpert. params)
1  153  1   1 143 1  1 80 1     > Active grid range: [iQTmin iQTmax iQTstep] [iymin ...] [iQmin ...]
1  153  1   1 143 1  1 80 1     > Full grid range (reference)
QT_______   y_______   Q______  > column header label (do not change)
```

**Key parameters to change:**

| Parameter | Description |
|-----------|-------------|
| `iBeam`   | 1=pp (LHC), -1=ppbar (Tevatron), 0=pN, -2=π⁻N |
| `ECM`     | √s in GeV (e.g. 8000 for LHC 8 TeV, 1960 for Tevatron) |
| `LTO`     | 0=full CSS+Y, -1=LO only, 3=Y-piece only, 1=ΔΣ integrated |
| `Type_V`  | W+, W-, Z0, H0, A0, AA, AG, ZZ, ZG |
| `ISTR`    | PDF set index (50=use LHAPDF via `lha_` prefix, 902=internal) |
| PDF name  | `lha_CT18NNLO`, `lha_CT14nnlo`, `lha_MSHT20nnlo_as118`, etc. |
| `iNONPERT`| 5=BLNY parametrization, 3=LY parametrization |
| `g1,g2,g3`| Non-perturbative Sudakov parameters |
| Grid range | Active rows/columns of the output grid to compute |

---

## SLURM Job Script for MSU HPCC

Save this as `run_legacy.sb` in your run directory:

```bash
#!/bin/bash --login
#SBATCH --job-name=legacy_W_LHC
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=8G
#SBATCH --partition=general-long

# Load compiler modules (must match what was used to compile)
module purge
module load GCC/12.3.0

# Set library paths
export PATH=/mnt/home/lopezels/InstallSources/LHAPDF/bin:$PATH
export LD_LIBRARY_PATH=/mnt/home/lopezels/InstallSources/LHAPDF/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/mnt/home/lopezels/InstallSources/HOPPET1/lib:$LD_LIBRARY_PATH

# OpenMP settings (Legacy uses OpenMP for grid integration)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_STACKSIZE=2G

# Move to the run directory
cd /mnt/home/lopezels/runs/legacy_W_LHC

# Run Legacy — output will be written to W_LHC_run.out
./main W_LHC_run

echo "Legacy finished with exit code $?"
```

Submit with:
```bash
sbatch run_legacy.sb
```

Monitor with:
```bash
squeue -u lopezels
tail -f slurm_<JOBID>.out
```

---

## Example: W⁺ Production at LHC 8 TeV

**Step 1.** Copy the code to a HPCC run directory:

```bash
mkdir -p /mnt/home/lopezels/runs/legacy_W_LHC8
cp -r /path/to/legacy_final_vesion/* /mnt/home/lopezels/runs/legacy_W_LHC8/
cd /mnt/home/lopezels/runs/legacy_W_LHC8
```

**Step 2.** Edit `legacy.in` for LHC 8 TeV W⁺:

```
1,0.d0,1,0                                       > iBeam=1 (pp), FRACT_N=0, KinCorr=1, ResumType=0
8000.0,0,0                                        > ECM=8000 GeV, LTO=0 (full CSS+Y), iFast=0
172., 80.358, 91.118, 125, 0                      > t,W,Z,h,photon masses
W+                                                > Boson type
50,2,1,2                                          > ISTR=50 (use LHAPDF), iPionPDF=2, i_RunMass=1
CSS                                               > Resummation scheme
5,99,1,1,1                                        > iNONPERT=5 (BLNY), IFLAG_C3=99
3,2,1,0,2                                         > Sudakov A=3, B=2, CFns=1, FSR=0, DelSigOrder=2
1,1,1                                             > iscale=1, muR=1, muF=1
lha_CT18NNLO                                      > Use CT18NNLO via LHAPDF6
0.06983d0,0.06984d0,0.15d0,1.3d0,4.75d0           > quark masses
0.5                                               > bMax
0                                                 > iProc=0 (all)
./inp/q_grid.inp                                  > Q grid
./inp/qt_grid.inp                                 > qT grid
./inp/y_grid.inp                                  > y grid
-
0.21d0, 0.68d0, -0.126d0, 1.6d0, 1               > g1,g2,g3,Q0,nG (BLNY)
1  153  1   1 143 1  1 80 1                       > Active grid
1  153  1   1 143 1  1 80 1                       > Full range
QT_______   y_______   Q______
```

**Step 3.** Compile (only needed once after code changes):

```bash
cd /mnt/home/lopezels/runs/legacy_W_LHC8

# Build HOPPET first (if not already done)
cd hoppet && ./configure --prefix=$(pwd) FC=gfortran && make -j4 && make install && cd ..

# Build Legacy
export PATH=/mnt/home/lopezels/InstallSources/LHAPDF/bin:$PATH
export LD_LIBRARY_PATH=/mnt/home/lopezels/InstallSources/LHAPDF/lib:$LD_LIBRARY_PATH
make clean && make
```

**Step 4.** Submit to SLURM:

```bash
sbatch run_legacy.sb
```

**Step 5.** Check output:

```bash
# After the job finishes:
head -100 W_LHC_run.out
```

The output file contains columns of `qT`, `y`, `Q`, and the differential cross section
`d²σ/dqT dy` in pb/GeV.

---

## Parallelism Notes

Legacy uses **OpenMP** to parallelize the integration over the (qT, y, Q) grid.
The number of threads is controlled by `OMP_NUM_THREADS`. On HPCC general nodes
you typically have 28 cores per node — matching your `.bashrc` setting of
`export OMP_NUM_THREADS=28`.

The `OMP_STACKSIZE=2G` is necessary because the CSS Bessel-transform integration
uses deep recursion that can overflow the default thread stack.

---

## Choosing a Partition on MSU HPCC

| Partition | Max walltime | Cores/node | Use when |
|-----------|-------------|------------|----------|
| `general` | 4 h | 28 | Quick test runs |
| `general-long` | 24 h | 28 | Full grid calculations |
| `bigmem` | 24 h | 28+ | Memory-intensive runs |

For a full (qT, y, Q) grid run with LTO=0, expect runtimes of **30 min – 4 h** depending
on grid density.

---

## Troubleshooting

**`lhapdf-config: command not found` during make**
Add LHAPDF to PATH before running make:
```bash
export PATH=/mnt/home/lopezels/InstallSources/LHAPDF/bin:$PATH
```

**`error while loading shared libraries: libLHAPDF.so.6`**
Add the LHAPDF lib directory at runtime:
```bash
export LD_LIBRARY_PATH=/mnt/home/lopezels/InstallSources/LHAPDF/lib:$LD_LIBRARY_PATH
```

**`fort.13` or `fort.2` files appear after a run**
These are diagnostic/warning files written by the code. `fort.13` contains quark mass
mismatch warnings (harmless — the code uses the value from `legacy.in`). `fort.2` contains
grid range errors. If you see entries in `fort.2`, check that your active grid indices
in `legacy.in` are within the bounds of the `.inp` files.

**OpenMP stack overflow / segfault**
Increase `OMP_STACKSIZE` in your SLURM script:
```bash
export OMP_STACKSIZE=4G
```

**PDF set not found**
Install the PDF set via LHAPDF:
```bash
/mnt/home/lopezels/InstallSources/LHAPDF/bin/lhapdf install CT18NNLO
```
