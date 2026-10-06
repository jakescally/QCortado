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

Auto compares one, two and four OpenMP threads per rank, with a modest
preference for two on large multi-node allocations. Other integer thread
counts are fallback choices when allocation divisibility, workflow limits or
source restart compatibility require them. It scores image/pool layouts by
node locality, alignment with node boundaries, point balance, and modest
replication/initialization overhead penalties. It favors pools occupying
roughly half to one node, rather than imposing a fixed ranks-per-pool target.
These weights are transparent heuristics informed by published benchmarks;
they are not a fitted performance model or timings measured on the user's cluster.

Known FFT dimensions and band counts from saved source SCF output provide
additional guidance for phonons. With pencil decomposition disabled, Auto
excludes layouts with more MPI ranks per pool than the smallest known FFT
third dimension, and penalizes approaching that limit. The phonon FFT grid
can differ from the source. When this constraint cannot be met, Auto asks for
a smaller allocation or a verified Manual layout. Missing dimensions fall
back to node locality. Auto does not predict memory requirements.

Phonon images distribute both q-points and irreps. **Known phonon irreps** in
Advanced settings optionally accepts the total selected irrep count reported
by QE initialization for the exact task, across all selected q-points. With
this count, Auto can choose more images than q-points, but never more than
available irreps: QE rejects images without an assigned representation.
Without it, Auto conservatively caps images at the estimated selected q-point
count. This is not a displacement-mode count or a count at just one q-point.
Changing the source, grid or selected representation ranges clears the supplied
count. Irrep costs differ, so balancing counts cannot guarantee balanced runtimes.

Pools use estimated electronic k-points; missing or conservative estimates use
one group unless preserving a recorded source SCF layout requires otherwise.
Compatible saved SCF ranks/pools per image take priority over performance
heuristics. Incompatible distributed restarts are flagged in Auto details.
Candidate counts must divide the MPI ranks. Auto preserves the requested CPU
budget and respects workflow rank limits. Pools spanning nodes or thread
counts above four trigger efficiency guidance in the standard tooltip.
For example, 20 nodes / 800 CPUs / 20 estimated electronic points without a
recorded source layout yields 400 ranks × 2 threads, with 20 pools for SCF.
A phonon task with 12 estimated q-points chooses 4 images × 5 pools per image,
20 MPI ranks per pool. Saved layouts and known FFT dimensions can change this.

The batch script exports `OMP_NUM_THREADS` from Slurm's CPUs per task,
propagates that setting through `SRUN_CPUS_PER_TASK`, uses block placement,
and defaults `OMP_PLACES=cores` / `OMP_PROC_BIND=close`, respecting preamble
overrides of those affinity settings. Profiles using `mpirun` should configure
their launcher arguments for thread-aware binding on that cluster.
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

Benchmark rationale: the [2023 QE developer study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10601483/)
includes phonon MPI/OpenMP and pool scaling, including a two-thread CPU layout;
the [ARCHER2 QE study](https://cug.org/proceedings/cug2022_proceedings/includes/files/pap104s2-file1.pdf)
shows that benefits from one, two and four threads depend strongly on the system.
[QE FFT scaling guidance](https://www.quantum-espresso.org/Doc/pw_user_guide/node18.html)
and [ARCHER2 launch guidance](https://docs.archer2.ac.uk/research-software/qe/)
inform FFT capacity and Slurm thread placement. None establishes a universally
optimal decomposition; compare representative timings on the target cluster.
