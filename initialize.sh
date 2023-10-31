#!/usr/bin/env bash

####################################################################################
# SCRIPT: initialize.sh
# DESCRIPTION:
# This messy script automates the setup for a Molecular Dynamics project. It performs the following tasks:
#   - Sets up the necessary folder structure.
#   - Installs and initializes Miniforge (Mambaforge) if necessary.
#   - Creates a Conda environment and installs required packages.
#   - Detects system information, including CPU, GPU, and CUDA versions.
#   - Installs additional libraries such as PyInteraph2 and MAVISp.
#   - Downloads starting structure files for specified proteins.
#   - Configures git and sets up a suitable .gitignore.
#   - test the MDAnalysis installation if the option -t
#
# PREREQUISITES:
#   - Assumes Miniforge (previously Mambaforge) is installed:
#       wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
#       bash Miniforge3-$(uname)-$(uname -m).sh
#       You should also activate the new environment using mamba init, although it can be done automatically
#   - Assumes CUDA and cuDNN are installed on the system (although currently not used).
#
# USAGE:
#   ./initialize.sh [PROTEIN_NAME]
#   PROTEIN_NAME: Name of the protein/system for which starting structure files should be downloaded.
#
# AUTHOR: Fabio Mazza
# DATE: 17/10/2023
# VERSION: 0.1
####################################################################################

## TODO:
#        fix genomicinfo.sh
#        support for multiple proteins or variants
#        make download_structures.sh more complete and flexible
#        personalize the folder structure based on the type of simulation(s) (REMD, GAMD etc.)
#        integrate with cookiecutter (!?)
#        change environment to inside project (with --prefix)?
#        make openmm and gamd packages/tests optional


# Initialize variables with default values
RUN_TESTS=false
INSTALL_OPENMM=false
ENVIRONMENT_NAME="md_project"  # Default environment name
ENV_NAME_SET=false # whether a new name has been passed to the script

# load author info
source .env
# set the yml files to use for the conda environment
CONDA_ENV_FILE="condaenv.yml"
CONDA_ENV_FILE_CUDA="condaenv_cuda.yml"

# Read the project folder path
PROJECT_FOLDER=$(pwd)


show_usage() {
    echo "Usage: $0 [-t] [-n ENVIRONMENT_NAME] [-c]"
    echo "  -t    Run tests after setting up the environment (optional)"
    echo "  -n    Specify the name of the conda environment (default: md_project)"
    echo "  -c    Use the CUDA-enabled conda environment file (condaenv_cuda.yml) instead of the default (condaenv.yml)"
    echo "  -o    Install and optionally test OpenMM and related packages"
    echo "  -h    Show this help message and exit"
}


# Function to read optional arguments
read_optional_args() {

    while getopts ":tn:coh" opt; do
        case ${opt} in
            t )
                RUN_TESTS=true
                ;;
            n )
                ENVIRONMENT_NAME="$OPTARG"
                ENV_NAME_SET=true
                ;;
            c )
                CONDA_ENV_FILE=$CONDA_ENV_FILE_CUDA
                ;;
            o )
                INSTALL_OPENMM=true
                ;;
            h )
                show_usage
                exit 0
                ;;
            \? )
                echo "Invalid Option: -$OPTARG" 1>&2
                show_usage
                exit 1
                ;;
            : )
                echo "Invalid Option: -$OPTARG requires an argument" 1>&2
                show_usage
                exit 1
                ;;
        esac
    done
    shift $((OPTIND -1))

    # Confirm with the user if they want to change the default environment name
    if [ "$ENVIRONMENT_NAME" = "md_project" ] && [ "$ENV_NAME_SET" = false ]; then
        echo -e "The default conda environment name is set to 'md_project'."
        read -p "Would you like to change it? (y/N): " change_name
        if [ "$change_name" = "y" ] || [ "$change_name" = "Y" ]; then
            echo -e "Please enter the new conda environment name:"
            read -r ENVIRONMENT_NAME
        fi
    fi


}



#####  DETECTION OF SYSTEM INFORMATION  #####

# Function to detect system information
detect_system_information() {
    echo
    echo "System Information:"
    echo -e "-------------------\n"

    # Kernel Version
    echo "Kernel Version:"
    uname -r
    echo

    # CPU Info
    echo -e "CPU Information:"
    echo
    echo -n "$(lscpu | grep "Model name" | awk '{$1=$1; print}') "
    avxlevel=$(./avxlevel.sh)
    numthreads=$(lscpu | grep "^CPU(s):" | awk '{print $2}')
    echo "(x86-64-v$avxlevel, $numthreads threads)"
    echo

    # GPU Info (Assumes NVIDIA GPU)
    if command -v nvidia-smi > /dev/null; then
        echo "GPU Information:"
        vram=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader)
        echo "$(nvidia-smi --query-gpu=name --format=csv,noheader) ($vram)"
    else
        echo "No GPU detected."
    fi
    echo

    # Check if nvcc is installed
    if command -v nvcc > /dev/null; then
        # CUDA Version
        echo "CUDA Version:"
        cuda_version=$(nvcc --version | grep "release" | awk '{print $6}' | cut -c2-)
        echo $cuda_version
        # CUDA Executable Location
        echo "CUDA Executable Location:"
        cuda_path=$(which nvcc)
        echo $cuda_path
    else
        echo "CUDA is not installed or nvcc is not in PATH"
    fi
    echo

    # Storage info
    DEVICE=$(df $PROJECT_FOLDER | tail -1 | awk '{print $1}')
    DEVICE_NAME=$(basename $DEVICE)
    DEVICE_NAME=${DEVICE_NAME%p*}
    DEVICE_TYPE=$(lsblk -d -o NAME,ROTA | grep $DEVICE_NAME | awk '{if ($2 == "0") print "SSD"; else print "HDD"}')
    DEVICE_FREE=$(df -h . | awk 'NR==2 {print $4 " free out of " $2}')

    echo -e "Storage Information:"
    echo "The storage of the project directory is a $DEVICE_TYPE with $DEVICE_FREE of space."

    FREE_SPACE_KiB=$(df --output=avail . | awk 'NR==2 {print $1}')
}



#####  Create folder structure  #####

create_folder_structure() {

    echo -e -n "\nCreating folder structure... "
    echo -e "-------------------\n"


    if (( FREE_SPACE_KiB < 300 * 1024 * 1024 )); then
        echo "Free space is less than 300GiB!"
        echo "It is recommended that you provide a new path with sufficient space for storing the project data. The data folder will be automatically created:"
        read -r NEW_PATH

        # Optionally, you might want to check if the provided path exists and has enough space
        while [[ ! -d "$NEW_PATH" ]]; do
            echo "The path does not exist. Please provide a valid path:"
            read -r NEW_PATH
        done

        DEVICE=$(df $NEW_PATH | tail -1 | awk '{print $1}')
        DEVICE_NAME=$(basename $DEVICE)
        DEVICE_NAME=${DEVICE_NAME%p*}
        DEVICE_TYPE=$(lsblk -d -o NAME,ROTA | grep $DEVICE_NAME | awk '{if ($2 == "0") print "SSD"; else print "HDD"}')
        DEVICE_FREE=$(df -h . | awk 'NR==2 {print $4 " free out of " $2}')

        echo "The storage of the project directory is an $DEVICE_TYPE with $DEVICE_FREE of space."

        mkdir -p "$NEW_PATH/data"

        # Create a symlink from the new path to your project data directory
        ln -s "$NEW_PATH" "$PROJECT_FOLDER"
        echo "Symlink created from $NEW_PATH to $PROJECT_FOLDER/data"
    else
        mkdir -p "data" # create data folder in the current
        echo -e "Data folder created in $PROJECT_FOLDER"
    fi

    mkdir -p "bin"
    mkdir -p "src"
    mkdir -p "scripts"
    mkdir -p "img"
    mkdir -p "papers"
    mkdir -p "data/00-structures"
    mkdir -p "data/01-minimization"
    mkdir -p "data/02-equilibration"
    mkdir -p "data/03-prod"
    mkdir -p "data/04-analysis"
    mkdir -p "data/forcefields"

    echo -e "Done\n\n"
}


#####  CREATE CONDA ENVIRONMENT  #####

create_conda_environment() {

    echo "Setup of the mamba environment"
    echo -e "-------------------\n"

    # update conda and install the faster mamba
    echo "Updating Conda and installing/upgrading Mamba..."
    conda update -n base -c conda-forge conda

    # Check if mamba is installed
    if ! command -v mamba &> /dev/null; then
        echo "Mamba is not installed. Installing now..."
        conda install mamba -c conda-forge --yes
    fi

    # Print the version of conda and Python
    echo "Mamba version:"
    mamba --version
    echo

    # Check if mamba is initialized
    # This checks for the conda initialize block in both .bashrc and .zshrc (common shells).
    # Adjust for other shells as needed.
    if ! grep -q "conda initialize" ~/.bashrc && ! grep -q "conda initialize" ~/.zshrc; then
        echo "Mamba is not initialized. Initializing now..."
        mamba init
        # Notify the user to restart the shell
        echo "Please restart your shell for the changes to take effect, then reexecute this script."
        exit 1
    fi

    # useless?
    CURRENT_SHELL=$(basename $SHELL)
    if [ "$CURRENT_SHELL" == "zsh" ]; then
        #source ~/.zshrc
        echo "skipping source ~/.zshrc"
    elif [ "$CURRENT_SHELL" == "bash" ]; then
        source ~/.bashrc
    # add more conditions for other shells if necessary
    else
        echo "Unsupported shell. Please restart your shell manually."
        exit 1
    fi


    # Paths to the new environment's folder, Python and pip executables
    #CONDA_PATH=$(conda info --envs | grep "^base " | awk '{print $2}') # 2 or 3?
    CONDA_PATH=$(conda info --envs | grep "^base " | awk '{print ($2 == "base" ? $3 : $2)}')
    # CONDA_PATH=$(conda config --show | grep "root_prefix" | awk '{print $2}') #alternative way

    # prefix path for the new environment
    ENV_PATH=$CONDA_PATH/envs/$ENVIRONMENT_NAME
    echo "Environment path: $ENV_PATH"

    # path for important executables, in case the environment cannot be activated from inside the script
    ENV_PYTHON_PATH=$CONDA_PATH/envs/$ENVIRONMENT_NAME/bin/python
    ENV_PIP_PATH=$(dirname $ENV_PYTHON_PATH)/pip
    echo "Environment python path: $ENV_PYTHON_PATH"

    # Set environment name in the yml files
    awk -v new_name="$ENVIRONMENT_NAME" 'NR==1 {sub(/name: .*/, "name: " new_name)} 1' $CONDA_ENV_FILE > temp.yml && mv temp.yml $CONDA_ENV_FILE
    awk -v new_name="$ENVIRONMENT_NAME" 'NR==1 {sub(/name: .*/, "name: " new_name)} 1' $CONDA_ENV_FILE_CUDA > temp.yml && mv temp.yml $CONDA_ENV_FILE_CUDA
    echo "Name of the environment set to $ENVIRONMENT_NAME."

    # Check if the INSTALL_OPENMM variable is set to false
    if [ "$INSTALL_OPENMM" == "false" ]; then
        # Use sed to delete the line containing "- openmm" from the environment files
        sed -i '/^- openmm/d' $CONDA_ENV_FILE
        sed -i '/^- openmm/d' $CONDA_ENV_FILE_CUDA
    fi

    echo -e "\nCreating the new conda environment..."
    export NVCC_ACCEPT_LICENSE=yes #automatic acceptance for CUDA stuff, if applicable

    set +e
    # Create a conda environment from the environment.yml file
    mamba env create -f "$CONDA_ENV_FILE" -p "$ENV_PATH" #--verbose
    set -e

    # update environment (useful for later script runs)
    echo "Updating the environment..."
    mamba env update -f "$CONDA_ENV_FILE" -p "$ENV_PATH" --prune

    echo "Activating the environment..."
    # Activate the environment (assuming mamba/conda is initialized in your shell)
    # This might not work from inside the script, we try it anyway
    #mamba activate $ENVIRONMENT_NAME
    # Alternative and old-ish way to activate the environment, it is slower and only works on Linux
    source $CONDA_PATH/etc/profile.d/conda.sh
    source $CONDA_PATH/etc/profile.d/mamba.sh
    mamba activate $ENVIRONMENT_NAME
    echo "Environment $ENVIRONMENT_NAME activated."

    # Set the project folder as a Conda environmental variable
    conda env config vars set PROJECT_DIR=$PROJECT_FOLDER # retrieve in python with project_dir = os.getenv("PROJECT_DIR")
    echo -e "Project folder set in the PROJECT_DIR Conda environmental variable.\n"

    # Enable the environment for notebooks
    $ENV_PYTHON_PATH -m ipykernel install --user --name $ENVIRONMENT_NAME

    echo
}


#####  INSTALL ADDITIONAL PACKAGES  #####

install_additional_packages() {

    # Extract the CUDA version in PATH
    # make sure that LD_LIBRARY_PATH is not set, since LD_LIBRARY_PATH can override the CUDA libraries
    # CUDA and CuDNN should be installed

    CUDA_VERSION=$(nvcc --version | grep "release" | awk '{print $6}' | cut -c2-3)

    # Determine which version of jax to install based on CUDA version
    if [[ $CUDA_VERSION == "12" ]]; then
        JAX_VERSION="jax[cuda12_pip]"
        echo "Installing JAX with CUDA $CUDA_VERSION backend..."
    elif [[ $CUDA_VERSION == "11" ]]; then
        JAX_VERSION="jax[cuda11_pip]"
        echo "Installing JAX with CUDA $CUDA_VERSION backend..."
    else
        JAX_VERSION="jax[cpu]"
        echo "Installing JAX without CUDA acceleration..."
    fi

    # Install the appropriate JAX version
    $ENV_PIP_PATH install --upgrade "$JAX_VERSION" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
    echo
    echo "Done."
    echo


    #conda install -c pytorch pytorch-cuda

    # FoldSeek local installation
    # static build with AVX2 (fastest)
    #wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz; tar xvfz mmseqs-linux-avx2.tar.gz; export PATH=$(pwd)/mmseqs/bin/:$PATH


    ## Check installed software

    if which vmd > /dev/null; then
        echo -e "\nVMD is installed."
    else
        echo "VMD is not installed or is not in PATH."
    fi

    if which pymol > /dev/null; then
        echo "pymol is installed."
    else
        echo "pymol is not installed or is not in PATH. Keep in mind an opensource and free version of PyMol can be found at github.com/schrodinger/pymol-open-source"
    fi

    if which gmx > /dev/null; then
        echo "GROMACS is installed."
    else
        echo -e "GROMACS is not installed or is not in PATH.\n You might want to add source /usr/local/gromacs/bin/GMXRC to your bashrc or zshrc.\n"
    fi
    echo



    ## Install additional optional libraries

    set +e
    echo -e "\nInstalling PyInteraph2..."
    cd ./bin
    git clone https://github.com/ELELAB/pyinteraph2.git
    cd pyinteraph2
    #pip install -r requirements.txt
    sed -i '/^\"MDAnalysis/d' setup.py # we disregard impositions on the minor version of MDAnalysis (and other libraries) by PyInteraph
    $ENV_PYTHON_PATH ./setup.py install
    echo
    cd ..

    echo "Installing MAVISp..."
    git clone https://github.com/ELELAB/MAVISp
    cd MAVISp
    $ENV_PYTHON_PATH install .
    echo
    cd ..

    echo "Downloading Prothon (structure ensemble comparison)..."
    git clone https://github.com/PlotkinLab/Prothon
    echo

    echo "Downloading ens-dRMS (structure ensemble comparison)..."
    git clone https://github.com/lazartomi/ens-dRMS
    echo

    # GAMD
    if [ "$INSTALL_OPENMM" = true ]; then
        git clone https://github.com/MiaoLab20/gamd-openmm
        cd gamd-openmm
        $ENV_PYTHON_PATH ./setup.py install
        cd ..
    fi

    git clone https://github.com/MiaoLab20/PyReweighting

    cd ..
    set -e
    echo
}


#####  DOWNLOAD STRUCTURES AND FORCEFIELDS  #####

download_structures() {
    echo -e "\nDownload starting structure and force field files..."
    echo -e "-------------------\n"

    cd data

    # Name of the protein/system
    echo -e "Please write the name of the gene you want to work on:"
    read NAME_GENE

    echo -e "Please write the name of the PDB of the wild-type protein, if it's available on RCSPDB:"
    read PDB_ID

    ../download_structures.sh "$NAME_GENE" "$PDB_ID"


    cd ./forcefields

    # Amber forcefields (most recent included in ambertools from the conda-forge channel, using the amber format)
    # see http://www.ks.uiuc.edu/Research/namd/2.9/ug/node13.html to use them with namd
    # Use parmed or CHARMM-GUI :( to convert them for gromacs use
    if [ ! -d "./Amber" ]; then
        mkdir "./Amber"
        if [ -d "$ENV_PATH/dat/leap" ]; then
            cp -r "$ENV_PATH/dat/leap/"* "./Amber/"
        else
            echo "The  directory $ENV_PATH/dat/leap does not exist. Ensure ambertools was correctly installed in your environment."
        fi
    else
        echo "The directory ./Amber already exists, skipping copy."
    fi

    # a99sb, best for disordered regions/proteins, for gromacs and amber https://doi.org/10.1101/2023.02.09.527891
    set +e
    git clone https://github.com/paulrobustelli/Force-Fields
    set -e

    cd ../..
}

#####  GIT  #####

configure_git() {

    echo -e "\n\nConfiguring git..."
    echo -e "-------------------\n"

    # List of trajectory extensions to ignore
    extensions=(
    ".xtc"
    ".trr"
    ".dcd"
    )
    folders=(
    "bin/"
    "/imgs/"
    "/papers/"
    "/data/"
    )

    # download and extend good gitignore, if not already present

    if [ ! -f "./.gitignore" ]; then
        curl https://www.toptal.com/developers/gitignore/api/python,julia > .gitignore
        # Append to .gitignore
        echo -e "\n## ADDITIONAL EXTENSIONS ##" >> .gitignore
        for ext in "${extensions[@]}"; do
            echo "$ext" >> .gitignore
        done
        echo -e "\n## PROJECT FOLDERS ##" >> .gitignore
        for fld in "${folders[@]}"; do
            echo "$fld" >> .gitignore
        done
        echo -e "\n.gitignore updated!\n"
        if [ $? -ne 0 ]; then
            echo -e "\nFailed to download .gitignore\n"
            exit 1
        fi
    else
        echo -e ".gitignore already exists, skipping download.\n"
    fi


    # initialize repo, if not already done
    if [ ! -d "./.git" ]; then
        git init
        git add .gitignore
        git add .
        git config user.name "$AUTHOR_NAME"
        git config user.email "$AUTHOR_EMAIL"
        git commit -m "Initial commit, folder structure and structure files added"
    else
        echo -e ".git already exists, skipping git init.\n"
    fi

    # Get the current Git branch
    current_branch=$(git rev-parse --abbrev-ref HEAD)

    # Check if the current branch is the default branch
    if [ "$current_branch" = "main" ] || [ "$current_branch" = "master" ]; then
    # Ask the user if they want to create a new branch
        read -p "Would you like to create a new branch? (y/N): " create_branch
        if [ "$create_branch" = "y" ] || [ "$create_branch" = "Y" ]; then
            echo -e "Please enter the name of the new branch:"
            read branch_name
            git checkout -b "$branch_name"
            echo -e "Switched to new branch: $branch_name\n"
        fi
    else
        echo -e "You are not on the default branch. Current branch: $current_branch\n"
    fi

    echo
}


#####   TEST PACKAGE INSTALLATION   #####

# if the -t argument is passed to initialize.sh, we performe all MDAnalysis concurrently
# only "fail" tests could pose problems, "xfail" can be disregarded. Plus, if any test fails,
# you should manually do a test without the --numprocesses option since it can sometimes cause problems

run_tests() {
    echo -e "\n\nRunning tests..."
    echo -e "-------------------\n"

    # Reactivate the environment (assuming mamba/conda is initialized in your shell)
    #mamba deactivate
    source $CONDA_PATH/etc/profile.d/conda.sh
    source $CONDA_PATH/etc/profile.d/mamba.sh
    set +e
    current_env=$(conda info --envs | grep '*' | awk '{print $1}')
    echo "Running tests in the Mamba/Conda Environment $current_env..."

    echo -e "Testing the MDAnalysis installation...\n"
    PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print tolower($1) substr($2, 1, 4)}')
    MDATPATH=$ENV_PATH/lib/$PYTHON_VERSION/site-packages/MDAnalysisTests/
    PYTEST_PATH="$ENV_PATH/bin/pytest"
    cd $MDATPATH
    $PYTEST_PATH  --disable-pytest-warnings --numprocesses $numthreads
    cd $PROJECT_FOLDER
    echo

    echo -e "\nJAX is using the following device:"

    # Use Python from the local environment to test that jax installed correctly with cuda support
    python3 ./test_python_packages.py

    # test openmm
    if [ "$INSTALL_OPENMM" == true ]; then
        echo -e "\nTesting OpenMM installation..."
        python3 -m openmm.testInstallation
    fi
    echo

}



main() {
    # Save the output of the script to a log file
    exec > >(tee -i initialize.log) 2>&1

    # Exit the script if any command fails
    set -e

    # make the necessary shell scripts executable
    chmod +x avxlevel.sh genomicinfo.sh download_structures.sh

    # Call functions
    read_optional_args "$@"
    # Print all optional arguments and their values
    echo "Optional Arguments and Their Values:"
    for arg in "${!OPTIONAL_ARGS[@]}"; do
        echo "$arg: ${OPTIONAL_ARGS[$arg]}"
    done

    detect_system_information
    create_folder_structure
    create_conda_environment
    install_additional_packages
    download_structures
    configure_git
    if [ "$RUN_TESTS" = true ]; then
        run_tests
    fi
}


# pass to main all command-line parameters and run it
main "$@"

