# Molecular Dynamics Project Initializer

Warning: still under active development.

This script automates the setup process for a Molecular Dynamics (MD) project, streamlining the initialization and configuration of necessary components.

## Features
- Sets up a sane default folder structure.
- Detects system information including CPU, GPU, and CUDA versions.
- Initializes Miniforge (Mambaforge) if not already set up.
- Creates a Conda environment and installs a good deal of useful packages for trajectory analysis.
- Searches and downloads starting structure files for specified proteins, from RCSPDB and/or the AlphaFold DB.
- Configures Git and sets up a complete .gitignore file.
- Runs installation tests (optional).

## Prerequisites
- Miniforge (previously Mambaforge) is recommended. You can download it from [here](https://github.com/conda-forge/miniforge#mambaforge).
- A CUDA-enabled GPU is recommended, even though the script currently does not set up GROMACS/NAMD/Amber.

## Usage
To use the script, navigate to your project directory and run:

```bash
./initialize.sh [options]
```

### Options

- `-t`: Run tests after setting up the environment.
- `-n ENVIRONMENT_NAME`: Specify the name of the conda environment (default: md_project).
- `-c`: Use the (condaenv_cuda.yml) conda environment file which includes CUDA libraries from the Nvidia conda channel.
- `-o`: Install OpenMM and OpenMM-related packages.
- `-h`: Show the help message and exit.


### Example

```bash
./initialize.sh -t -n my_md_project
```

This example runs the initialization script, performs tests after setting up the environment, and specifies "my_md_project" as the conda environment name.

## Folder Structure

The script will create a series of directories to organize your MD project:

- `bin`: For binary files.
- `src`: For additional packages with source code.
- `scripts`: For additional scripts.
- `img`: For image files and plots.
- `papers`: For related papers and documentation.
- `data`: For data files, further organized into subdirectories for different stages of the MD project, excluded from git versioning.

