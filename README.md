# Modeling contact binaries
## III. Properties of a population of massive, close binaries

These files accompany the article by:
* Matthias Fabry,
* Pablo Marchant,
* Norbert Langer,
* Hugues Sana,

submitted to Astronomy and Astrophysics for publication.
With them, you can reproduce the results of the paper.

The files fall, broadly, in two categories.
1. Creating and running the mesa models, and
2. analysing and plotting the results.

### Category 1
* `zams_models/` contains the non-rotating ZAMS models that the grid starts from. Copy the contents of this folder into your MESA installation folder in `data/star_data/zams_models`
* `make_grid.py` uses `makefolders.py` to build the directory structure in which to run the MESA models that use energy transfer/no energy transfer.
It copies from `templates/grid`, and uses `run_both`.
If you put your new `grid` folder at the same level, you should not have to change any paths in any file.
`run_grid.slurm` is there for your convenience to run the models easily on a cluster. Issuing `sbatch -a 1-n run_grid.slurm` will start `n` parallel jobs, running the models of the grid. You will have to change pathnames here for your system.
* `make_single_rot_grid.py` uses `makefolders.py` to build the directories in which the single-rotating models will be run.
`run_single_rot_grid.slurm` is used to run the grid calculations on a cluster (mind path differences for your system).

To make a grid using the `grid` mesa template, issue
`python make_grid.py grid_folder templates/grid`
which will create the `grid_folder/` directory containing the appropriate folder structure for the parameter variations you specified.
To run a single model, `cd` into one, _e.g._ `cd grid_folder/28.0/1.26/0.775/ET/` and issue  
`../../../../../../execs/mk`  
`../../../../../../execs/rnrt`  
(cf. the commands in `run_both`). The reason for this pedantic executable calling is to avoid having lots of duplicate files in the directory tree, to save file count when using cluster resources.
Same for the `inlists` and other files in the `common_files/` directory.

* the full grid of models (`grid.tar.gz`) is available on Zenodo (doi:10.5281/zenodo.14289889), if you got this repository from there, it should be included. This file is ~30GB and highly compressed, so it might take a while to uncompress depending on your system resources (you might want to consider `pigz`). The uncompressed size is **~250GB**.
* `cases/` contains the models needed to reproduce figures 1-3, they are copied from the `grid/` and included here so you do not have to download the full grid tar file.
* `appendices/` contains the models that were computed for Appendices B and C, and come with their full history information included.
* `data_tables/` contains the data required to use the tidal deformation ocrrections, computed in Fabry et al. 2022 (A&A,Vol.661,A123)

### Category 2
#### Helper files
* `mesa_data.py`, for loading `history.data`
* `mesa_population.py`, for holding a bunch of histories in one object
* `pop_synth.py`, functions for doing population synthesis
* `observational_data.py`, parameters of observed contact systems

#### Plotting files
* `make_plots_calculations.py`, main file where calculations and most of the plot making is issued
* The `npzs/` directory contains files that story the properties of the population for fast access during plotting. 
Recomputing the population would take an appreciable amount of time loading all the models.
* `completion.py`, makes the large completion matrix figures
* `plot_example_models.py`, makes the plots of the example models of Class I, II and IV, as well those featured in the Appendices
* `RLOF_at_zams_figures.py`, makes the figure of models that overflow at the zams.


