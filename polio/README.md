## Polio lesson

The polio lesson exemplifies a more sophisticated workflow for a real scientific problem.

----------------------------

### Run-levels

The workflow illustrated in this lesson makes use of "run-levels", as controlled by an **R** variable `run_level`.
The choice of run-level affects various algorithmic parameters in such a way that higher values of the run-level result in more extensive (and expensive), but more thorough calculations.
At run-level 1, the codes are quickly completed: this is useful for debugging.
At run-level 2, more computational resources are mobilized: this is useful for checking that the computations will scale.
Run-level 3 is intended to produce meaningful results.

In the codes herein, one changes the run-level by setting the environment variable `RUNLEVEL` to 1, 2, or 3.
Alternatively, one can edit the codes to set the run-level explicitly.

The slides and notes are generated using run-level 3: the directory `results` contains the archived computational results.
The results of running the codes in `main.R` at run-levels 1 and 2 are contained in the directories `runlevel1`, `runlevel2`, respectively.


----------------------------

### Manifest

- main.Rnw, main.R: main lesson
- algorithmic-parameters-exercise.html
- convergence-exercise.html
- demography-exercise.html
- initial-values-exercise.html
- starting-values-exercise.html
- params.csv: parameter estimates
- polio_wisconsin.csv: data
- subdirectories:
  - results: archive files for `main.Rnw` (RUNLEVEL=3)
  - runlevel1: results of running `main.R` with RUNLEVEL=1
  - runlevel2: results of running `main.R` with RUNLEVEL=2
  - initial-values-exercise: worked solution
  - starting-values-exercise: worked solution

----------------------------

### License

![CC-BY-NC](https://i.creativecommons.org/l/by-nc/4.0/88x31.png)

This work is licensed under the [Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0)](http://creativecommons.org/licenses/by-nc/4.0/).
Under its terms, you are free to:

- Share — copy and redistribute the material in any medium or format
- Adapt — remix, transform, and build upon the material

The licensor cannot revoke these freedoms as long as you follow the license terms.
Under the following terms:

- Attribution — You must give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.
- NonCommercial — You may not use the material for commercial purposes.
- No additional restrictions — You may not apply legal terms or technological measures that legally restrict others from doing anything the license permits.

----------------------------
