# Quantum ESPRESSO point parallelization

In the SCF MPI settings or the HPC run settings for SCF, bands and DOS,
set **K-point pools** to split the electronic calculation using QE's `-nk`
option. Phonon settings additionally offer **Phonon images**, using
`ph.x -nimage` to distribute q-points and irreducible representations.
In Manual mode both counts default to 1; explicitly configured point layouts
retain their manual settings.

CPU cluster runs offer **Auto** within the existing HPC allocation controls.
Enter **Nodes** and **Total allocated CPUs** (the total across all nodes).
Auto fills MPI tasks, CPUs per task, placement, phonon images and electronic
pools, and shows the resulting breakdown beneath the allocation fields.
Partition, walltime, memory and advanced cluster settings remain in the same
settings block. Switch to Manual to edit the decomposition directly; GPU
allocations retain manual controls.
Auto disables pencil decomposition and locks its checkbox; Manual permits
enabling it explicitly for compatible cluster installations.

Auto evaluates MPI/OpenMP splits together with image and pool layouts.
Allocations of at least 128 CPUs prefer four threads per rank when there are
at least 32 CPUs per node, two with at least 16 CPUs per node, and one
otherwise. It can choose other integer thread counts to preserve the CPU
budget, fit node placement and avoid oversized MPI groups. It favors phonon
images using the estimated irreducible q-count, then electronic pools using
the estimated k-count, aiming for roughly eight MPI ranks per pool. Eight is a
starting heuristic, not a QE requirement or a guaranteed optimum. Candidate
counts must divide the ranks, with a preference for balanced point workloads.
When possible phonons preserve the saved source SCF ranks/pools per image;
incompatible distributed restarts are identified in the Auto details tooltip.
Missing or conservative point estimates use one group for that level.
Workflow rank limits also constrain candidate layouts. Auto always preserves
the requested total CPU budget or reports an invalid allocation. Layouts with
more than eight threads per rank or more than 64 CPUs per electronic pool
flag potentially excessive allocation in the Auto details tooltip. These
thresholds are conservative heuristics, not measured scaling limits.
For example, 20 nodes / 800 CPUs / 20 irreducible electronic points yields
160 ranks × 5 threads, with 20 pools of 8 ranks each. Actual choices also
depend on phonon point counts and saved source SCF MPI layouts.
The batch script exports `OMP_NUM_THREADS` from Slurm's CPUs per task and
defaults `OMP_PLACES=cores` / `OMP_PROC_BIND=close`, respecting preamble
overrides of those affinity settings. The MPI launcher must allow each rank
to use its allocated CPUs; profiles using `mpirun` should configure their
launcher arguments for thread-aware binding on that cluster.
Evenly divisible allocations set tasks per node; uneven allocations let Slurm
place ranks rather than retaining a conflicting profile default.

The SCF and DOS grid controls preview the expected irreducible k-point count.
Phonons preview irreducible q-points and the source SCF k-mesh. The pool controls
offer balanced counts that divide both the estimated point count and MPI tasks
per image. Band paths show their explicit sampled count, including endpoints
and disconnected segments, rather than a uniform-mesh symmetry estimate.
Guidance is available through the standard information tooltips.

Mesh previews use the matching calculation cell (including the inferred QE
Bravais basis), grid offsets and spglib symmetry reduction. They are estimates:
QE's tolerance and magnetic symmetries can change the final result. For
noncollinear/SOC settings the preview gives a conservative full-mesh bound
instead of claiming a magnetic symmetry reduction. Phonon k-point symmetry
also changes with q; the source SCF count is guidance, not a strict pool cap.

The number of MPI tasks must be divisible by pools × images. CPUs per task
are OpenMP threads, not additional MPI ranks. Use no more pools than the
number of irreducible k-points, preferably a divisor of that number. Avoid
excessive ranks per pool relative to FFT dimensions. Keeping a pool within
a node can reduce communication; the best layout depends on the system and
cluster, so compare timings before scaling up.

For example, a 64-task phonon allocation with 8 images and 2 pools uses
8 tasks per image and 4 tasks per pool. QE's documented compatible source
SCF layout is 8 tasks with 2 pools. SCF k-points share and update the same
charge density through QE; they must not be submitted as independent SCF jobs.

On clusters, the preview and submitted batch script use the same decomposition.
Both configured launchers (`srun` and `mpirun`) are supported. The allocation
size stays as requested. Image collection runs with one image and the original
tasks per image / pools, sets `recover = .true.`, and completes before `q2r.x`
and `matdyn.x`. Collection failure stops the pipeline. Local MPI phonons also
perform this recovery pass automatically.

Use a shared filesystem for the SCF save data and all phonon image scratch.
Keep the image/pool layout unchanged when recovering an interrupted run, and
retain its scratch files. Distributed-format QE restart files require the same
process decomposition that wrote them. The wizard records the requested
decomposition with saved SCF and phonon calculation parameters.

References: [QE parallelization levels](https://www.quantum-espresso.org/Doc/user_guide/node20.html),
[phonon parallelism and collection](https://www.quantum-espresso.org/Doc/ph_user_guide/node17.html),
[parallel I/O](https://www.quantum-espresso.org/Doc/user_guide/node21.html).
