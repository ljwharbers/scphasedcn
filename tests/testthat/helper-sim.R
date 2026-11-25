# helper-sim.R: shared simulated data sets for the test suite (built once per session).
# sim_small: tumour regime with multiplets; sim_flat: one flat cell type, no events, no multiplets.
# seed 1: the data these tests were written on (see helper-caller.R; in sim_small the allele-error draws shifted).
sim_small = simulate_scn(n_cells = 160, seed = 1)
sim_flat = simulate_scn(n_cells = 120, seed = 1, multiplet_rate = 0, events = NULL, cell_types = c(normal = 1))
flat_cells = function(sim) sim$scn$cells$cell[!sim$truth$cells$multiplet & sim$scn$cells$cell_type == "normal"]
clean_cells = function(sim) sim$scn$cells$cell[!sim$truth$cells$multiplet]
