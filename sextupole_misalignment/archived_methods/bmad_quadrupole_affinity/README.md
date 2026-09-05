# Historical Bmad/Tao quadrupole-affinity reference

Archived on 2026-09-05 from `quadrupole_affinity/`. This is an explicitly labeled
Bmad/Tao reference calculation, retained for historical method comparison and
Bmad-versus-SciBmad validation. The primary latest-lattice SciBmad study remains
in [`../../quadrupole_affinity/`](../../quadrupole_affinity/README.md).

## Method and preserved files

The historical generator screens quadrupole K1 interventions using local
beta/phase changes, then calculates finite-difference K2/offset response
dictionaries. It includes responses to the target offset and nominal nuisance
responses to the other sextupole offsets. This predates the maintained SciBmad
direct-trajectory response dictionary.

- `generate_bmad_affinity_responses.py`: historical PyTao response generator.
- `results/responses/`: saved Bmad response bundles and generation metadata.
- `results/affinity/`: historical scores, figures, and comparison summaries.
- `validate_nominal_responses.py`: compares saved nominal Bmad responses with
  the archived SciBmad/GTPSA mixed-derivative coefficients in
  [`../response_map/`](../response_map/README.md).

All result files and their original metadata were moved unchanged. Embedded
absolute paths record where those files were originally generated.

## Validate the saved cross-code comparison

From `CESR Project/`, using Python with NumPy:

```bash
python -B sextupole_misalignment/archived_methods/bmad_quadrupole_affinity/validate_nominal_responses.py --output-dir=/tmp/cesr-bmad-affinity-validation
```

This reads existing Bmad and SciBmad data; it does not launch Tao or generate a
new lattice scan. The migration check reproduced the saved comparison tables.

The generator is retained for historical reproduction and explicit cross-code
diagnosis only. Any new Bmad/PyTao execution must use the project's documented
`Ubuntu-Bmad` WSL2 distribution and `bmad` Conda environment. Its project-root
resolution now accounts for the archive directory; the latest SciBmad affinity
generator and analysis utilities remain in the maintained study.
