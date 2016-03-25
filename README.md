# nse

A modified version of the nse package ( https://r-forge.r-project.org/projects/nse/ ) that adds bootstrappers for xts objects based on the existing cpp bootstrapper via hashing the dates and passing them, then unhashing them.

Note that the wrapper to apply functions on the bootstrapped samples does not work with most PerfA metrics as they are bugged when they encounter multiple samples with the same date tag , which bootstrapping often causes.

Hence, the wrappers should be taken as a proof of concept, and that alone.
