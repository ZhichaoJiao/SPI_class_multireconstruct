#!/bin/bash

size=100
n_pattern=35000
lambda=1e-10
z_det=0.5
pix_len=600e-6
step=0.1
r_min=3
r_max=32
n_gamma=200
initial_mode_1=1    #1:input volume, 2:random orientation
volume_1_path="/ssd1/Data/SPARTA/sparta_pattern_3/pattern_3_volume/volume_8isz_af2_resize_100_2.h5"
initial_mode_2=1
volume_2_path="/ssd1/Data/SPARTA/sparta_pattern_3/pattern_3_volume/volume_8k9g_af2_resize_100_2.h5"
initial_mode_3=1
volume_3_path="/ssd1/Data/SPARTA/sparta_pattern_3/pattern_3_volume/volume_8it1_af2_resize_100_2.h5"
pattern_path="/ssd1/Data/SPARTA/sparta_pattern_3/pattern_mix_1/pattern_resize_100_2/pattern_"
output_path="/ssd1/Data/SPARTA/sparta_pattern_3/pattern_mix_1/output_2"

cc_dif=0.02
iter=100


###########以下不用修改##########
cp $0 ${output_path}/parameter.txt

mkdir "${output_path}/iteration_1"

#orientation_1
/usr/lib64/mpich/bin/mpirun \
    -np 8 \
    ./orientation_one_iter \
    --size=$size \
    --pattern_path=$pattern_path \
    --volume_path=$volume_1_path \
    --output_path="${output_path}/iteration_1" \
    --n_pattern=$n_pattern \
    --lambda=$lambda \
    --z_det=$z_det \
    --pix_len=$pix_len \
    --step=$step \
    --r_min=$r_min --r_max=$r_max \
    --n_gamma=$n_gamma \
    --fine_search=0 \
    --index_volume=1 \
    --initial_mode=$initial_mode_1

#orientation_2
/usr/lib64/mpich/bin/mpirun \
    -np 8 \
    ./orientation_one_iter \
    --size=$size \
    --pattern_path=$pattern_path \
    --volume_path=$volume_2_path \
    --output_path="${output_path}/iteration_1" \
    --n_pattern=$n_pattern \
    --lambda=$lambda \
    --z_det=$z_det \
    --pix_len=$pix_len \
    --step=$step \
    --r_min=$r_min --r_max=$r_max \
    --n_gamma=$n_gamma \
    --fine_search=0 \
    --index_volume=2 \
    --initial_mode=$initial_mode_2

#orientation_3
/usr/lib64/mpich/bin/mpirun \
    -np 8 \
    ./orientation_one_iter \
    --size=$size \
    --pattern_path=$pattern_path \
    --volume_path=$volume_3_path \
    --output_path="${output_path}/iteration_1" \
    --n_pattern=$n_pattern \
    --lambda=$lambda \
    --z_det=$z_det \
    --pix_len=$pix_len \
    --step=$step \
    --r_min=$r_min --r_max=$r_max \
    --n_gamma=$n_gamma \
    --fine_search=0 \
    --index_volume=3 \
    --initial_mode=$initial_mode_3

for i in $(seq 1 ${iter})
do
    output_path_i="${output_path}/iteration_${i}"
    output_path_i_next="${output_path}/iteration_$((i+1))"
    mkdir $output_path_i_next
    #compare_cc cc_list_volume_1:所有衍射图与volume_1切片计算cc结果. pattern_list_volume_1:用于merge volume_1的所有衍射图

    /home/jiao/anaconda3/envs/orientation/bin/python \
    ./compare_cc_3.py \
    --infile_cc_1="${output_path_i}/cc_list_volume_1.h5" \
    --infile_cc_2="${output_path_i}/cc_list_volume_2.h5" \
    --infile_cc_3="${output_path_i}/cc_list_volume_3.h5" \
    --cc_dif=$cc_dif \
    --outfile_cc_1="${output_path_i}/pattern_list_volume_1.h5" \
    --outfile_cc_2="${output_path_i}/pattern_list_volume_2.h5" \
    --outfile_cc_3="${output_path_i}/pattern_list_volume_3.h5" 

    #merge_volume_1
    ./merge_volume \
    --size=$size \
    --lambda=$lambda \
    --z_det=$z_det \
    --pix_len=$pix_len \
    --pattern_path=$pattern_path \
    --output_file="${output_path_i}/merge_volume_1.h5" \
    --angle_path="${output_path_i}/pattern_list_volume_1.h5" 

    #merge_volume_2
    ./merge_volume \
    --size=$size \
    --lambda=$lambda \
    --z_det=$z_det \
    --pix_len=$pix_len \
    --pattern_path=$pattern_path \
    --output_file="${output_path_i}/merge_volume_2.h5" \
    --angle_path="${output_path_i}/pattern_list_volume_2.h5" 

    #merge_volume_3
    ./merge_volume \
    --size=$size \
    --lambda=$lambda \
    --z_det=$z_det \
    --pix_len=$pix_len \
    --pattern_path=$pattern_path \
    --output_file="${output_path_i}/merge_volume_3.h5" \
    --angle_path="${output_path_i}/pattern_list_volume_3.h5" 

    #orientation_1
    /usr/lib64/mpich/bin/mpirun \
        -np 8 \
        ./orientation_one_iter \
        --size=$size \
        --pattern_path=$pattern_path \
        --volume_path="${output_path_i}/merge_volume_1.h5" \
        --output_path=$output_path_i_next \
        --n_pattern=$n_pattern \
        --lambda=$lambda \
        --z_det=$z_det \
        --pix_len=$pix_len \
        --step=$step \
        --r_min=$r_min --r_max=$r_max \
        --n_gamma=$n_gamma \
        --fine_search=0 \
        --index_volume=1 \
        --initial_mode=1

    #orientation_2
    /usr/lib64/mpich/bin/mpirun \
        -np 8 \
        ./orientation_one_iter \
        --size=$size \
        --pattern_path=$pattern_path \
        --volume_path="${output_path_i}/merge_volume_2.h5" \
        --output_path=$output_path_i_next \
        --n_pattern=$n_pattern \
        --lambda=$lambda \
        --z_det=$z_det \
        --pix_len=$pix_len \
        --step=$step \
        --r_min=$r_min --r_max=$r_max \
        --n_gamma=$n_gamma \
        --fine_search=0 \
        --index_volume=2 \
        --initial_mode=1

    #orientation_3
    /usr/lib64/mpich/bin/mpirun \
        -np 8 \
        ./orientation_one_iter \
        --size=$size \
        --pattern_path=$pattern_path \
        --volume_path="${output_path_i}/merge_volume_3.h5" \
        --output_path=$output_path_i_next \
        --n_pattern=$n_pattern \
        --lambda=$lambda \
        --z_det=$z_det \
        --pix_len=$pix_len \
        --step=$step \
        --r_min=$r_min --r_max=$r_max \
        --n_gamma=$n_gamma \
        --fine_search=0 \
        --index_volume=3 \
        --initial_mode=1
done








