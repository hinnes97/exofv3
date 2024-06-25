### Basic setup of GCM
- Recursively clone this repo, `git clone --recurse-submodules -j8 git@github.com:hinnes97/exofv3.git`
- Move to the `tools` directory and run `build_postprocessing.bash`, i.e. `cd tools && ./build_postprocessing.bash`. Running this script requires GNU autotools and the necessary environment variables populated so that autotools can find the necessary compilers, netcdf and openmpi libraries required for postprocessing. Currently these are populated using `module load` commands at the beginning of the script, which are set up to be used on the Oxford servers. To run on a non-Oxford system, these will need to be deleted or replaced.
- Ensure there is a mkmf template file appropriate for your system in `tools/mkmf_templates`. This will tell the mkmf makefile generator where to find the compilers/netcdf/mpi libraries used in the compilation of the main code. An example Oxford template can be found in this directory
- Change directory to the default experiment `exofv3/exp/default`
- Create a file named `user_config` with the same format as the `user_config_example` file provided. This will be read in by `run.bash` and sets some of the main options of the code, like the output directory and the grid resolution.
- Try running the code, using `./run.bash`. The runscript should take care of compilation, execution and postprocessing if configured correctly.
