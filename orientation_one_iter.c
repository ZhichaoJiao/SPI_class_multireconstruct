/*
只迭代一轮,计算并输出所有衍射图的cc_max,不进行merge
使用orientation_main 20231120 版本
  */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <hdf5.h>
#include <sys/io.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include "function.h"
#define pi 3.141592653590

int size;
int n_gamma;
int n_r;

int main(int argc, char *argv[])
{
    // 初始化mpi
    MPI_Init(NULL, NULL);
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // 声明变量
    // 可调参数
    int fine_search = 0; // 是否进行精细搜索, 精细搜索的步长为step_fine = step/fine_search
    int N_pattern = 0;   // 衍射图总数
    double lambda = 0;   // X射线波长，单位 m
    double z_det = 0;    // 探测器与样品距离，单位 m
    double pix_len = 0;  // 像素点尺寸， 单位 m
    double step = 0;     // 三维采样间隔,单位 弧度
    int r_min = 0;
    int r_max = 0;
    int index_volume = 0; // volume编号
    char pattern_path[400];
    char initial_volume_path[400];
    char output[400];
    int initial_mode=1; //初始模式,1.输入volume,2.随机取向

    // 从命令行输入参数
    const struct option longopts[] = {
        {"n_pattern", 1, NULL, 0},
        {"lambda", 1, NULL, 2},
        {"z_det", 1, NULL, 3},
        {"pix_len", 1, NULL, 4},
        {"step", 1, NULL, 5},
        {"r_min", 1, NULL, 6},
        {"r_max", 1, NULL, 7},
        {"n_gamma", 1, NULL, 8},
        {"volume_path", 1, NULL, 11},
        {"pattern_path", 1, NULL, 12},
        {"output_path", 1, NULL, 13},
        {"size", 1, NULL, 15},
        {"fine_search", 1, NULL, 16},
        {"index_volume", 1, NULL, 17},
        {"initial_mode", 1, NULL, 18},
        {0, 0, NULL, 0}};

    /*short options*/
    /*从命令中读取参数信息，如果有输入，那么将输入赋给参数，如果没有输入，那么使用参数变量声明时候的初始值*/
    int c;
    char *rval;
    while ((c = getopt_long(argc, argv, "a", longopts, NULL)) != -1)
    {

        switch (c)
        {
        case 0:
            N_pattern = strtod(optarg, &rval);
            break;

        case 2:
            lambda = strtod(optarg, &rval);
            break;

        case 3:
            z_det = strtod(optarg, &rval);
            break;

        case 4:
            pix_len = strtod(optarg, &rval);
            break;

        case 5:
            step = strtod(optarg, &rval);
            break;

        case 6:
            r_min = strtol(optarg, &rval, 10);
            break;

        case 7:
            r_max = strtol(optarg, &rval, 10);
            break;

        case 8:
            n_gamma = strtol(optarg, &rval, 10);
            break;

        case 11:
            snprintf(initial_volume_path, 400, "%s", optarg);
            break;

        case 12:
            snprintf(pattern_path, 400, "%s", optarg);
            break;

        case 13:
            snprintf(output, 400, "%s", optarg);
            break;

        case 15:
            size = strtol(optarg, &rval, 10);
            break;

        case 16:
            fine_search = strtol(optarg, &rval, 10);
            break;

        case 17:
            index_volume = strtol(optarg, &rval, 10);
            break;

        case 18:
            initial_mode = strtol(optarg, &rval, 10);
            break;
        }
    }

    // 输出参数
    if (0 == world_rank)
    {
        printf("\n");
        printf("size=%d\n", size);
        printf("pattern_path : %s\n", pattern_path);
        printf("initial_volume_path : %s\n", initial_volume_path);
        printf("output_path : %s\n", output);
        printf("N_pattern = %d\n", N_pattern);
        printf("lambda = %e\n", lambda);
        printf("z_det = %f\n", z_det);
        printf("pix_len = %e\n", pix_len);
        printf("step = %f\n", step);
        printf("r_min = %d\n", r_min);
        printf("r_max = %d\n", r_max);
        printf("n_gamma = %d\n", n_gamma);
        printf("fine_search = %d\n", fine_search);
        printf("initial_mode = %d\n", initial_mode);
        printf("\n");
        char paramater_file_path[400];
        FILE *paramater_file;
        snprintf(paramater_file_path, 400, "%s/paramater.txt", output);
        paramater_file = fopen(paramater_file_path, "w");
        fprintf(paramater_file, "initial_mode:(1)import_volume;(2)random_angle\n");
        fprintf(paramater_file, "size = %d\n", size);
        fprintf(paramater_file, "pattern_path : %s\n", pattern_path);
        fprintf(paramater_file, "initial_volume_path : %s\n", initial_volume_path);
        fprintf(paramater_file, "output_path : %s\n", output);
        fprintf(paramater_file, "N_pattern = %d\n", N_pattern);
        fprintf(paramater_file, "lambda = %e\n", lambda);
        fprintf(paramater_file, "z_det = %f\n", z_det);
        fprintf(paramater_file, "pix_len = %e\n", pix_len);
        fprintf(paramater_file, "step = %f\n", step);
        fprintf(paramater_file, "r_min = %d\n", r_min);
        fprintf(paramater_file, "r_max = %d\n", r_max);
        fprintf(paramater_file, "n_gamma = %d\n", n_gamma);
        fprintf(paramater_file, "fine_search = %d\n", fine_search);
        fclose(paramater_file);
    }

    // 其他变量
    n_r = r_max - r_min;
    double rot_angle[3] = {0, 0, 0};          // 衍射图旋转角,rot_angle={alpha, beta, gamma}.以入射X射线方向为z轴负方向，建立右手坐标系，rot_angle表示这个坐标系相对于样品颗粒坐标系的旋转角度，三个角度依次对应于zxz顺规下的欧拉角。
    double best_angle_matrix_i[N_pattern][5]; // 用于保存每个进程中,所有衍射图找到的最佳角度.第一列为衍射图序号,每行的五个参数分别是{Index_pattern, alpha, beta, gamma, cc}
    memset(best_angle_matrix_i, 0, sizeof(best_angle_matrix_i));
    double(*pattern_cart)[size] = (double(*)[size])malloc(sizeof(double) * size * size); // 直角坐标系下二维衍射图强度二维数组
    memset(pattern_cart, 0, sizeof(double) * size * size);
    float(*pattern_cart_float)[size] = (float(*)[size])malloc(sizeof(float) * size * size); // float类型的衍射图,用来读取hdf5文件
    memset(pattern_cart_float, 0, sizeof(float) * size * size);
    double reference_cart[size][size];
    memset(reference_cart, 0, sizeof(reference_cart));
    double reference_polar[n_r][n_gamma];
    memset(reference_polar, 0, sizeof(reference_polar));
    float *volume_1d = (float *)malloc(sizeof(float) * size * size * size); // 三维衍射强度一维数组
    memset(volume_1d, 0, sizeof(float) * size * size * size);
    double *pixels = (double *)malloc(sizeof(double) * size * size * 4); // 二维衍射图的坐标数据，已经变化到Ewald sphere上面，一维double数组，pixels = {x1, y1, z1, omega1, ..., }
    memset(pixels, 0, sizeof(double) * size * size * 4);

    int n_angle_max = 0.1 * 0.1 * 1400 / (step * step); // 根据步长估算全空间三维采样点数量(最大值),步长为0.1弧度时,采样点~1400个
    int N_reference = 0;                                // 切片数量
    double rot_angle_matrix[n_angle_max][3];            // 记录空间取向的矩阵，每一行代表一个取向，三列依次代表zxz顺规下的三个欧拉角[alpha, beta, gamma]，单位是弧度。
    memset(rot_angle_matrix, 0, sizeof(rot_angle_matrix));

    // 生成采样角度
    MakeAllAngle(rot_angle_matrix, step, &N_reference);
    int batch = N_reference / world_size + 1; // 每个进程分配的切片数量
    double rot_angle_matrix_i[batch][3];      // 每个进程分配的切片对应的取向
    memset(rot_angle_matrix_i, 0, sizeof(double) * batch * 3);
    double(*reference_polar_matrix_i)[n_r][n_gamma] = (double(*)[n_r][n_gamma])malloc(sizeof(double) * n_r * n_gamma * batch); // 用来保存为每个进程分配的所有切片
    memset(reference_polar_matrix_i, 0, sizeof(double) * n_r * n_gamma * batch);
    double(*best_angle_matrix)[N_pattern][5]; // 用于保存所有进程汇总后,所有衍射图的角度,只在0进程中分配内存和使用,大小为world_size * N_pattern * 5. 第n行即为第n张衍射图.第一列为衍射图序号,每行的五个参数分别是{Index_pattern, alpha, beta, gamma, cc}

    // 为每个进程分配取向
    for (int i = 0; i < batch; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int t = world_rank * batch;
            rot_angle_matrix_i[i][j] = rot_angle_matrix[t + i][j];
        }
    }

    // 计算出各个像素点对应于Ewald球上的坐标
    GenPixels(pixels, lambda, z_det, pix_len, size);

    // 生成初始Volume, 单进程
    if (0 == world_rank)
    {
        printf("N_reference = %d, batch = %d\n", N_reference, batch);
        printf("make initial volume...\n");

        double real_angle[N_pattern][4];
        memset(real_angle, 0, sizeof(double) * N_pattern * 4);
        double(*volume_3d)[size][size] = (double(*)[size][size])malloc(sizeof(double) * size * size * size); // 三维衍射强度三维数组
        memset(volume_3d, 0, sizeof(double) * size * size * size);
        float(*volume_3d_float)[size][size] = (float(*)[size][size])malloc(sizeof(float) * size * size * size); // 三维衍射强度三维数组
        memset(volume_3d_float, 0, sizeof(float) * size * size * size);
        double *weight = (double *)malloc(sizeof(double) * size * size * size); // 每个三维像素点的总权重
        memset(weight, 0, sizeof(double) * size * size * size);
        best_angle_matrix = (double(*)[N_pattern][5])malloc(sizeof(double) * N_pattern * 5 * world_size);
        memset(best_angle_matrix, 0, sizeof(double) * N_pattern * 5 * world_size);
        if (1 == initial_mode)
        {
            // 输入initial volume
            hid_t file_volume_ini, dataset, dataspace;
            herr_t status;
            file_volume_ini = H5Fopen(initial_volume_path, H5F_ACC_RDONLY, H5P_DEFAULT);
            dataset = H5Dopen(file_volume_ini, "volume_pow", H5P_DEFAULT);
            status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, volume_3d_float);
            H5Dclose(dataset);
            status = H5Fclose(file_volume_ini);

            Volume3dFloat2Double(volume_3d_float, volume_3d);
        }

        else if (2 == initial_mode)
        {
            // 按照随机取向Merge
            srand((unsigned)time(NULL)); // 随机数种子
            for (int index_pattern = 0; index_pattern < N_pattern; index_pattern++)
            {
                char pattern_path_i[400] = {0};
                memset(pattern_cart, 0, sizeof(double) * size * size);
                memset(pattern_cart_float, 0, sizeof(float) * size * size);

                // 生成随机取向
                double random_number_1 = rand() % 10000 / 10000.0; // 随机数范围[0,1],精确到万分位小数,注意除数必须是浮点数类型的10000.0,不能写成整数类型的10000
                double random_number_2 = rand() % 10000 / 10000.0;
                double random_number_3 = rand() % 10000 / 10000.0;
                rot_angle[0] = 2 * pi * random_number_1;
                rot_angle[1] = pi * random_number_2;
                rot_angle[2] = 2 * pi * random_number_3;

                // 读取衍射图数据
                snprintf(pattern_path_i, 400, "%s%d.h5", pattern_path, index_pattern);
                hid_t file_pattern, dataset, dataspace;
                herr_t status;
                file_pattern = H5Fopen(pattern_path_i, H5F_ACC_RDONLY, H5P_DEFAULT);
                dataset = H5Dopen(file_pattern, "data", H5P_DEFAULT);
                status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pattern_cart_float);
                H5Dclose(dataset);
                status = H5Fclose(file_pattern);

                PatternCartFloat2Double(pattern_cart_float, pattern_cart);

                PatternMerge(rot_angle, pattern_cart, volume_1d, weight, size, pixels);
            }
            // 按照权重weight将volume_1d转化成volume_3d
            Volume_1d2Volume_3d(volume_3d, volume_1d, weight, size);
        }

        // 将volume_3d转化为volume_1d. 注意volume_1d转化成volume_3d的过程中加入了权重weight,所以两个转化过程不可逆,也不可省略
        Volume_3d2Volume_1d(volume_3d, volume_1d);

        free(volume_3d);
        free(volume_3d_float);
        free(weight);
    }

    // 将进程0中计算得到的volume_1d广播到所有进程
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(volume_1d, size * size * size, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // 每个进程独立生成各自的切片
    if (0 == world_rank)
        printf("Make reference...\n");

    for (int index_reference = 0; index_reference < batch; index_reference++)
    {
        memset(reference_cart, 0, sizeof(double) * size * size);
        memset(reference_polar, 0, sizeof(double) * n_r * n_gamma);

        rot_angle[0] = rot_angle_matrix_i[index_reference][0];
        rot_angle[1] = rot_angle_matrix_i[index_reference][1];
        rot_angle[2] = rot_angle_matrix_i[index_reference][2];
        ReferenceGen(rot_angle, reference_cart, volume_1d, size, pixels);                             // 切片函数
        Cart2Polar(reference_cart, reference_polar_matrix_i[index_reference], r_min, r_max, n_gamma); // 转换成极坐标
    }

    // for循环,每张衍射图寻找在各自进程中的最佳取向
    char pattern_path_i[400] = {0};
    double cc_max = 0;
    double(*pattern_polar)[n_gamma] = (double(*)[n_gamma])malloc(sizeof(double) * n_r * n_gamma);
    memset(pattern_polar, 0, sizeof(double) * n_r * n_gamma);
    double best_angle[3] = {0, 0, 0};
    double(*dif_angle_cc_i)[4] = (double(*)[4])malloc(sizeof(double) * batch * 4); // 衍射图在这个进程中角度的cc系数分布
    memset(dif_angle_cc_i, 0, sizeof(double) * batch * 4);
    double best_angle_cc[n_gamma];
    memset(best_angle_cc, 0, sizeof(double) * n_gamma);
    double(*dif_angle_cc)[4] = (double(*)[4])malloc(sizeof(double) * (batch * world_size) * 4); // 储存前三张衍射图在所有角度的cc系数分布
    memset(dif_angle_cc, 0, sizeof(double) * (batch * world_size) * 4);

    for (int index_pattern = 0; index_pattern < N_pattern; index_pattern++)
    {
        cc_max = 0;
        memset(pattern_cart, 0, sizeof(double) * size * size);
        memset(pattern_cart_float, 0, sizeof(float) * size * size);
        memset(pattern_polar, 0, sizeof(double) * n_r * n_gamma);
        memset(dif_angle_cc_i, 0, sizeof(double) * batch * 4);
        memset(best_angle_cc, 0, sizeof(double) * n_gamma);
        memset(dif_angle_cc, 0, sizeof(double) * (batch * world_size) * 4);

        // 读取衍射图数据
        snprintf(pattern_path_i, 400, "%s%d.h5", pattern_path, index_pattern);
        hid_t file_pattern_i, dataset_i;
        herr_t status_i;
        file_pattern_i = H5Fopen(pattern_path_i, H5F_ACC_RDONLY, H5P_DEFAULT);
        dataset_i = H5Dopen(file_pattern_i, "data", H5P_DEFAULT);
        status_i = H5Dread(dataset_i, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pattern_cart_float);
        H5Dclose(dataset_i);
        status_i = H5Fclose(file_pattern_i);

        PatternCartFloat2Double(pattern_cart_float, pattern_cart);

        Cart2Polar(pattern_cart, pattern_polar, r_min, r_max, n_gamma);
        FindBestAngle(pattern_polar, reference_polar_matrix_i, rot_angle_matrix_i, best_angle, batch, r_min, r_max, n_gamma, dif_angle_cc_i, best_angle_cc, &cc_max);

        best_angle_matrix_i[index_pattern][0] = index_pattern;
        best_angle_matrix_i[index_pattern][1] = best_angle[0];
        best_angle_matrix_i[index_pattern][2] = best_angle[1];
        best_angle_matrix_i[index_pattern][3] = best_angle[2];
        best_angle_matrix_i[index_pattern][4] = cc_max;

        // 输出第一张衍射图的cc分布
        if (1 == index_pattern || 200001 == index_pattern)
        {
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Gather(dif_angle_cc_i, batch * 4, MPI_DOUBLE, dif_angle_cc, batch * 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            if (0 == world_rank)
            {
                double cc_max = 0;
                double best_angle[4] = {0};
                for (int index_cc = 0; index_cc < N_reference; index_cc++)
                {
                    if (dif_angle_cc[index_cc][3] > cc_max)
                    {
                        cc_max = dif_angle_cc[index_cc][3];
                        best_angle[0] = dif_angle_cc[index_cc][0];
                        best_angle[1] = dif_angle_cc[index_cc][1];
                        best_angle[2] = dif_angle_cc[index_cc][2];
                        best_angle[3] = dif_angle_cc[index_cc][3];
                    }
                }
                hid_t file_out_cc;
                herr_t status;
                hid_t dataspace, dataset;
                hsize_t dim_1[2] = {batch * world_size, 4}, dim_2[2] = {1, 4};
                char file_out_cc_path[400] = {0};
                snprintf(file_out_cc_path, 400, "%s/cc_pattern_%d_volume_%d.h5", output, index_pattern, index_volume);
                file_out_cc = H5Fcreate(file_out_cc_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
                dataspace = H5Screate_simple(2, dim_1, NULL);
                dataset = H5Dcreate(file_out_cc, "dif_angle_cc", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dif_angle_cc);
                H5Dclose(dataset);
                H5Sclose(dataspace);
                dataspace = H5Screate_simple(2, dim_2, NULL);
                dataset = H5Dcreate(file_out_cc, "best_angle", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, best_angle);
                H5Dclose(dataset);
                H5Sclose(dataspace);
                status = H5Fclose(file_out_cc);
            }
        }

        if (0 == world_rank)
        {
            if (index_pattern < 10)
            {
                printf("    pattern_%d/%d cc_max=%f\n", index_pattern, N_pattern, cc_max);
            }
            else if (index_pattern % 100 == 0)
            {
                printf("\r    pattern_%d/%d cc_max=%f", index_pattern, N_pattern, cc_max);
                fflush(stdout);
            }
            if (index_pattern == N_pattern - 1)
            {
                printf("\n    pattern_%d/%d cc_max=%f\n", index_pattern, N_pattern, cc_max);
            }
        }
    }

    // 收集所有衍射图的最佳角度到0号进程
    if (0 == world_rank)
        printf("Waiting for other process...\n");
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(best_angle_matrix_i, N_pattern * 5, MPI_DOUBLE, best_angle_matrix, N_pattern * 5, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // bug2
    if (0 == world_rank)
    {
        // 判断MPI_Gather后的数据格式是否正确
        for (int a = 0; a < world_size; a++)
        {
            for (int b = 0; b < N_pattern; b++)
            {
                if (best_angle_matrix[a][b][0] != b)
                {
                    printf("\nError: Bug 2\nbest_angle_matrix[%d][%d][%d]=%f\n", a, b, 0, best_angle_matrix[a][b][0]);
                    return 1;
                }
            }
        }

        // 判断MPI_Gather前后数据是否一致
        for (int b = 0; b < N_pattern; b++)
        {
            for (int c = 0; c < 5; c++)
            {
                if (best_angle_matrix[0][b][c] != best_angle_matrix_i[b][c])
                {
                    printf("\nError: Bug 3\nbest_angle_matrix[0][%d][%d]=%f, best_angle_matrix_i[%d][%d]=%f, not equal\n", b, c, best_angle_matrix[0][b][c], b, c, best_angle_matrix_i[b][c]);
                    return 1;
                }
            }
        }
    }

    // 整合所有进程的最佳取向,找到每张衍射图的最佳取向,将最终结果保存在0进程中的best_angle_matrix_i中
    if (0 == world_rank)
    {
        printf("Caculate globle cc_max in process 0...\n");
        for (int index_pattern = 0; index_pattern < N_pattern; index_pattern++)
        {
            double cc_max = 0;
            for (int index_rank = 0; index_rank < world_size; index_rank++)
            {
                if (best_angle_matrix[index_rank][index_pattern][4] > cc_max)
                {
                    cc_max = best_angle_matrix[index_rank][index_pattern][4];
                    best_angle_matrix_i[index_pattern][0] = best_angle_matrix[index_rank][index_pattern][0];
                    best_angle_matrix_i[index_pattern][1] = best_angle_matrix[index_rank][index_pattern][1];
                    best_angle_matrix_i[index_pattern][2] = best_angle_matrix[index_rank][index_pattern][2];
                    best_angle_matrix_i[index_pattern][3] = best_angle_matrix[index_rank][index_pattern][3];
                    best_angle_matrix_i[index_pattern][4] = best_angle_matrix[index_rank][index_pattern][4];
                }
            }
        }
    }

    if (fine_search)
    {
        // 精细搜索:对每个衍射图,在最佳取向附近两个步长范围进行精细搜索,搜索步长缩小十倍
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(best_angle_matrix_i, N_pattern * 5, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        int batch_pattern = N_pattern / world_size + 1;
        double step_fine = step / fine_search;
        int n_fine_search = (fine_search * 2) * (fine_search * 2);
        double rot_angle_matrix_fine[n_fine_search][3]; // 精细搜索过程中,每个进程内每张衍射图最佳取向附近的所有精细搜索角度
        memset(rot_angle_matrix_fine, 0, sizeof(double) * n_fine_search * 3);
        double best_angle_matrix_fine_i[batch_pattern][5]; // 精细搜索过程中,每个进程内所有衍射图的最佳取向与cc系数
        memset(best_angle_matrix_fine_i, 0, sizeof(double) * batch_pattern * 5);
        double(*reference_polar_matrix_fine_i)[n_r][n_gamma] = (double(*)[n_r][n_gamma])malloc(sizeof(double) * n_r * n_gamma * n_fine_search); // 精细搜索过程中,每个进程内的所有切片
        memset(reference_polar_matrix_fine_i, 0, sizeof(double) * n_r * n_gamma * n_fine_search);
        double(*dif_angle_cc_fine_i)[4] = (double(*)[4])malloc(sizeof(double) * n_fine_search * 4); // 衍射图在这个进程中角度的cc系数分布
        memset(dif_angle_cc_fine_i, 0, sizeof(double) * n_fine_search * 4);

        if (0 == world_rank)
            printf("Fine searching: step_fine = %f, N_reference = %d\n", step_fine, n_fine_search);

        for (int index_pattern = world_rank * batch_pattern; index_pattern < min((world_rank + 1) * batch_pattern, N_pattern); index_pattern++)
        {
            best_angle[0] = best_angle_matrix_i[index_pattern][1];
            best_angle[1] = best_angle_matrix_i[index_pattern][2];
            best_angle[2] = best_angle_matrix_i[index_pattern][3];

            MakeFineAngle(best_angle, rot_angle_matrix_fine, step, step_fine, fine_search);

            for (int index_reference = 0; index_reference < n_fine_search; index_reference++)
            {
                memset(reference_cart, 0, sizeof(reference_cart));
                memset(reference_polar, 0, sizeof(reference_polar));

                rot_angle[0] = rot_angle_matrix_fine[index_reference][0];
                rot_angle[1] = rot_angle_matrix_fine[index_reference][1];
                rot_angle[2] = rot_angle_matrix_fine[index_reference][2];
                ReferenceGen(rot_angle, reference_cart, volume_1d, size, pixels);                                  // 切片函数
                Cart2Polar(reference_cart, reference_polar_matrix_fine_i[index_reference], r_min, r_max, n_gamma); // 转换成极坐标
            }

            // 读取衍射图数据
            snprintf(pattern_path_i, 400, "%s%d.h5", pattern_path, index_pattern);
            hid_t file_pattern_i, dataset_i;
            herr_t status_i;
            file_pattern_i = H5Fopen(pattern_path_i, H5F_ACC_RDONLY, H5P_DEFAULT);
            dataset_i = H5Dopen(file_pattern_i, "data", H5P_DEFAULT);
            status_i = H5Dread(dataset_i, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pattern_cart_float);
            H5Dclose(dataset_i);
            status_i = H5Fclose(file_pattern_i);

            PatternCartFloat2Double(pattern_cart_float, pattern_cart);
            Cart2Polar(pattern_cart, pattern_polar, r_min, r_max, n_gamma);
            FindBestAngle(pattern_polar, reference_polar_matrix_fine_i, rot_angle_matrix_fine, best_angle, n_fine_search, r_min, r_max, n_gamma, dif_angle_cc_fine_i, best_angle_cc, &cc_max);
            int index = index_pattern - world_rank * batch_pattern;
            best_angle_matrix_fine_i[index][0] = index_pattern;
            best_angle_matrix_fine_i[index][1] = best_angle[0];
            best_angle_matrix_fine_i[index][2] = best_angle[1];
            best_angle_matrix_fine_i[index][3] = best_angle[2];
            best_angle_matrix_fine_i[index][4] = cc_max;

            if (0 == world_rank)
            {
                if (index_pattern < 10)
                {
                    printf("    pattern_%d/%d cc_max=%f\n", index_pattern, batch_pattern, cc_max);
                }
                else if (index_pattern % 100 == 0)
                {
                    printf("\r    pattern_%d/%d cc_max=%f", index_pattern, batch_pattern, cc_max);
                    fflush(stdout);
                }
                if (index_pattern == batch_pattern - 1)
                {
                    printf("\n    pattern_%d/%d cc_max=%f\n", index_pattern, batch_pattern, cc_max);
                }
            }
        }

        // 收集所有衍射图的最佳角度到0号进程
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gather(best_angle_matrix_fine_i, batch_pattern * 5, MPI_DOUBLE, best_angle_matrix_i, batch_pattern * 5, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        free(reference_polar_matrix_fine_i);
        free(dif_angle_cc_fine_i);
    }

    // 输出所有衍射图的最佳角度
    if (0 == world_rank)
    {
        hid_t file_out_volume;
        herr_t status;
        hid_t dataspace, dataset;
        char file_out_volume_path[400] = {0};
        snprintf(file_out_volume_path, 400, "%s/cc_list_volume_%d.h5", output, index_volume);
        file_out_volume = H5Fcreate(file_out_volume_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        // best_angle_matrix
        hsize_t dim_4[2] = {N_pattern, 5};
        dataspace = H5Screate_simple(2, dim_4, NULL);
        dataset = H5Dcreate(file_out_volume, "angle", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, best_angle_matrix_i);
        H5Dclose(dataset);
        H5Sclose(dataspace);
        status = H5Fclose(file_out_volume);
    }

    free(pattern_polar);
    free(dif_angle_cc_i);
    free(dif_angle_cc);

    free(pattern_cart);
    free(volume_1d);
    free(pixels);
    free(reference_polar_matrix_i);
    free(pattern_cart_float);

    if (0 == world_rank)
    {
        printf("\n----All Finish !!!----\n");
        free(best_angle_matrix);
    }

    MPI_Finalize();
    return 0;
}