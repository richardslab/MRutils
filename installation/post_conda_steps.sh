#source me

ENV=vitaminD_MR

conda activate ${ENV}

## run post-conda steps
echo RUNNING post-conda steps

R --no-save < installation/post_conda_steps.R

