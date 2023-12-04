echo 'stack.sh = sbatch GMX files'

read -p 'first ndx: ' first_ndx
read -p 'last ndx: ' last_ndx

cd ./[...] # EDIT - workdir 
for ndxout in $( seq $first_ndx $last_ndx )
do
    let 'ndxin = ndxout - 1'
    SHFILE=[...]$ndxout.sh # EDIT - rename [...] (ex. cas9_$ndxout.sh)
    cat > $SHFILE << stack_eof
#!/bin/bash
#SBATCH --job-name cas9_$ndxout
#SBATCH -N1 --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=09:55:00
#SBATCH --gres=gpu:1
#SBATCH --error=stack.err
#SBATCH --output=stack.txt

module load gromacs/2020.4
export OMP_NUM_THREADS=16

gmx grompp -f [...].mdp -c [...].gro -t [...]$ndxin.cpt -p [...].top -n [...].ndx -o [...]$ndxout.tpr # EDIT [...]
gmx mdrun -s [...]$ndxout.tpr -v -ntmpi 1 -deffnm [...]$ndxout # EDIT [...]
gmx trjconv -f [...]$ndxout.trr -s [...]$ndxout.tpr -n [...].ndx -pbc cluster -o [...]$ndxout.xtc << zero_awk # EDIT [...]
0
0
zero_awk
stack_eof

    stack_id=$(squeue -u glattanzi | grep "cas9_$ndxin" | awk '{ print $3 }')
    if [ -z "$stack_id" ]
    then
        sbatch $SHFILE
    else
        sbatch --dependency=afterok:$stack_id $SHFILE
    fi
done
