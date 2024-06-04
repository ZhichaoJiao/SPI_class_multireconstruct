#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include <sys/io.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <libgen.h>

#define pi 3.141592653590
int size;

void GenPixels(double *pixels, double lambda, double z_det, double pix_len, int size);
void PatternMerge(double rot_angle[3], double pattern_cart[size][size], float *volume_1d, double *weight, int size, double *pixels);
void MakeRotMatrixEuler(double *rot_angle, double rot_matrix[3][3]);
void PatternCartFloat2Double(float pattern_cart_float[size][size], double pattern_cart_double[size][size]);
void Volume_1d2Volume_3d(double volume_3d[size][size][size], float *volume_1d, double *weight, int size);
void Volume3dDouble2FLoat(double volume_3d_double[size][size][size], float volume_3d_float[size][size][size]);

int main(int argc, char *argv[])
{
	size = 100;
	double lambda = 1e-10;
	double z_det = 1;
	double pix_len = 1e-6;
    char angle_path[400];
	char pattern_path[400];
	char output_file[400];

   // 从命令行输入参数

    const struct option longopts[] = {
        {"lambda", 1, NULL, 2},
        {"z_det", 1, NULL, 3},
        {"pix_len", 1, NULL, 4},
        {"pattern_path", 1, NULL, 5},
        {"output_file", 1, NULL, 6},
        {"angle_path", 1, NULL, 7},
        {"size", 1, NULL, 8},
        {0, 0, NULL, 0}};

	printf("Merge volume begin...\n");

    /*short options*/
    /*从命令中读取参数信息，如果有输入，那么将输入赋给参数，如果没有输入，那么使用参数变量声明时候的初始值*/
    int c;
    char *rval;
    while ((c = getopt_long(argc, argv, "a", longopts, NULL)) != -1)
    {

        switch (c)
        {

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
            snprintf(pattern_path, 400, "%s", optarg);
            break;

        case 6:
            snprintf(output_file, 400, "%s", optarg);
            break;

        case 7:
            snprintf(angle_path, 400, "%s", optarg);
            break;

        case 8:
            size = strtol(optarg, &rval, 10);
            break;
        }
    }

	double pattern_cart[size][size];
	float pattern_cart_float[size][size];
	double rot_angle[3];
	float *volume_1d = (float *)malloc(sizeof(float) * size * size * size); // 三维衍射强度一维数组
	memset(volume_1d, 0, sizeof(float) * size * size * size);
	double *weight = (double *)malloc(sizeof(double) * size * size * size); // 每个三维像素点的总权重
	memset(weight, 0, sizeof(double) * size * size * size);
	double *pixels = (double *)malloc(sizeof(double) * size * size * 4); // 二维衍射图的坐标数据，已经变化到Ewald sphere上面，一维double数组，pixels = {x1, y1, z1, omega1, ..., }
	memset(pixels, 0, sizeof(double) * size * size * 4);
	double(*volume_3d)[size][size] = (double(*)[size][size])malloc(sizeof(double) * size * size * size); // 三维衍射强度三维数组
	memset(volume_3d, 0, sizeof(double) * size * size * size);
	float(*volume_3d_float)[size][size] = (float(*)[size][size])malloc(sizeof(float) * size * size * size); // 三维衍射强度三维数组
	memset(volume_3d_float, 0, sizeof(float) * size * size * size);

	GenPixels(pixels, lambda, z_det, pix_len, size);

	hid_t file_real_angle, dataset, dataspace;
	herr_t status;
	hsize_t dims[2];
	file_real_angle = H5Fopen(angle_path, H5F_ACC_RDONLY, H5P_DEFAULT);
    dataset = H5Dopen2(file_real_angle, "angle", H5P_DEFAULT);
    dataspace = H5Dget_space(dataset);
    status = H5Sget_simple_extent_dims(dataspace, dims, NULL);

	int N_pattern=dims[0];
	double real_angle[N_pattern][5];
	memset(real_angle,0 ,sizeof(double)*N_pattern*5);
	printf("N_pattern = %d\n",N_pattern);

	status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, real_angle);
	H5Dclose(dataset);
	status = H5Fclose(file_real_angle);

	char *path=dirname(strdup(output_file));
	if (access(path, W_OK) == 0) {
        printf("Output path is right\n");
		printf("%s\n",output_file);
    } else {
		printf("Error: Output path is wrong\n");
		printf("%s\n",output_file);
		return 0;
    }
	
	for (int index_angle = 0; index_angle < N_pattern; index_angle++)
	{	

		int index_pattern;
		index_pattern=real_angle[index_angle][0];

		char pattern_path_i[400] = {0};
		memset(pattern_cart, 0, sizeof(double) * size * size);
		memset(pattern_cart_float, 0, sizeof(float) * size * size);

		// 输入正确取向
		rot_angle[0] = real_angle[index_angle][1];
		rot_angle[1] = real_angle[index_angle][2];
		rot_angle[2] = real_angle[index_angle][3];

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

	// 输出Merge_volume
	Volume3dDouble2FLoat(volume_3d, volume_3d_float);
	hid_t file_merge_volume;
	char file_out_volume_path[400] = {0};
	file_merge_volume = H5Fcreate(output_file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	hsize_t dim_3[3] = {size, size, size};
	dataspace = H5Screate_simple(3, dim_3, NULL);
	dataset = H5Dcreate(file_merge_volume, "volume_pow", H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, volume_3d_float);
	H5Dclose(dataset);
	H5Sclose(dataspace);
	status = H5Fclose(file_merge_volume);

	printf("\nFinish ! ! !\n");

	return 0;
}

void GenPixels(double *pixels, double lambda, double z_det, double pix_len, int size)
{
	// 根据模拟衍射的各项参数，构造出每个像素点对应的坐标数组Pixels = {x1, y1, z1, omega1, x2, ..., }，以倒空间坐标原点为原点，向下为x轴正方形，向右为y轴正方向，建立右手系
	// 计算结果:探测器上x方向的最大分辨率,也是Volume中x方向最大的分辨率,Volume中y,z方向的最大分辨率与x相同.

	double x_plane, y_plane, z_plane, x_sphere, y_sphere, z_sphere; // plane代表探测器平面上像素点的坐标，sphere代表投影到Ewald球后曲面上像素点的坐标。
	double length;													// 以正空间样品位置为原点的坐标, 探测器与样品之间的距离length.
	double pixels0;													// 储存(x0^2+y0^2)^(1/2)，用于坐标的归一化
	double omega, cos_2theta;										// 立体角omega，衍射角2theta
	long i, j, t;
	double center, dsize;
	dsize = size;			  // 将size转化为double类型
	center = (dsize - 1) / 2; // 中心位置

	// 计算出pixels0，物理含义为x=0，y=center探测器平面上的点，变换到Ewald球后得到的x坐标。
	i = 0;
	t = 0;
	x_plane = (i - center) * pix_len;
	y_plane = 0;
	z_plane = -z_det; // 为了满足建立的右手坐标系，z值需要取负
	length = sqrt(x_plane * x_plane + y_plane * y_plane + z_plane * z_plane);
	x_sphere = x_plane / (length * lambda);
	pixels0 = fabs(x_sphere);

	// 计算pixels[]，坐标参数缩放到（-center，center）范围
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			t = i * size + j;
			x_plane = (i - center) * pix_len;
			y_plane = (j - center) * pix_len;
			z_plane = -z_det;
			length = sqrt(x_plane * x_plane + y_plane * y_plane + z_plane * z_plane);
			x_sphere = x_plane / (length * lambda);
			y_sphere = y_plane / (length * lambda);
			z_sphere = z_plane / (length * lambda);
			cos_2theta = z_det / length;
			omega = (pix_len * pix_len) * cos_2theta / (length * length);
			pixels[t * 4] = center * x_sphere / pixels0;
			pixels[t * 4 + 1] = center * y_sphere / pixels0;
			pixels[t * 4 + 2] = center * (z_sphere + 1 / lambda) / pixels0;
			pixels[t * 4 + 3] = omega;
		}
	}
}

void PatternMerge(double rot_angle[3], double pattern_cart[size][size], float *volume_1d, double *weight, int size, double *pixels)
{
	// 根据衍射图的旋转角，将一张衍射图的强度插值到三维，得到Merge三维衍射强度的一维数组volume_1d和新的权重数组weight
	//  size表示每个维度（x/y）的像素点个数
	long x, y, z, N_pix = size * size;
	double tx, ty, tz, fx, fy, fz, cx, cy, cz, w, f;
	double rot_pix[3], rot_matrix[3][3] = {{0}};
	double *slice = (double *)malloc(sizeof(double) * size * size); // 直角坐标系下二维切片的一维数组
	memset(slice, 0, sizeof(double) * size * size);
	double center, dsize;
	dsize = size;			  // 将size转化为double类型
	center = (dsize - 1) / 2; // 中心位置

	// 构造旋转矩阵
	MakeRotMatrixEuler(rot_angle, rot_matrix);

	// 根据二维衍射图数据，将二维衍射强度pattern_cart变成一个一维数组slice
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			long t = i * size + j;
			slice[t] = pattern_cart[i][j];
		}
	}

	// 根据旋转矩阵，求出每个像素点旋转后的三维空间坐标，记录在rot_pix[3]中
	for (long t = 0; t < N_pix; t++)
	{
		for (int i = 0; i < 3; i++)
		{
			rot_pix[i] = 0.;
			for (int j = 0; j < 3; j++)
			{
				rot_pix[i] += rot_matrix[i][j] * pixels[t * 4 + j];
			}
			rot_pix[i] += center;
		}

		tx = rot_pix[0];
		ty = rot_pix[1];
		tz = rot_pix[2];

		x = tx;
		y = ty;
		z = tz;

		if (x < 1 || x > size - 2 || y < 1 || y > size - 2 || z < 1 || z > size - 2)
			continue;

		fx = tx - x;
		fy = ty - y;
		fz = tz - z;
		cx = 1. - fx;
		cy = 1. - fy;
		cz = 1. - fz;
		if (0 != slice[t])
		{
			slice[t] /= pixels[t * 4 + 3];
			w = slice[t];

			f = cx * cy * cz;
			weight[x * size * size + y * size + z] += f;
			volume_1d[x * size * size + y * size + z] += f * w;

			f = cx * cy * fz;
			weight[x * size * size + y * size + ((z + 1) % size)] += f;
			volume_1d[x * size * size + y * size + ((z + 1) % size)] += f * w;

			f = cx * fy * cz;
			weight[x * size * size + ((y + 1) % size) * size + z] += f;
			volume_1d[x * size * size + ((y + 1) % size) * size + z] += f * w;

			f = cx * fy * fz;
			weight[x * size * size + ((y + 1) % size) * size + ((z + 1) % size)] += f;
			volume_1d[x * size * size + ((y + 1) % size) * size + ((z + 1) % size)] += f * w;

			f = fx * cy * cz;
			weight[((x + 1) % size) * size * size + y * size + z] += f;
			volume_1d[((x + 1) % size) * size * size + y * size + z] += f * w;

			f = fx * cy * fz;
			weight[((x + 1) % size) * size * size + y * size + ((z + 1) % size)] += f;
			volume_1d[((x + 1) % size) * size * size + y * size + ((z + 1) % size)] += f * w;

			f = fx * fy * cz;
			weight[((x + 1) % size) * size * size + ((y + 1) % size) * size + z] += f;
			volume_1d[((x + 1) % size) * size * size + ((y + 1) % size) * size + z] += f * w;

			f = fx * fy * fz;
			weight[((x + 1) % size) * size * size + ((y + 1) % size) * size + ((z + 1) % size)] += f;
			volume_1d[((x + 1) % size) * size * size + ((y + 1) % size) * size + ((z + 1) % size)] += f * w;
		}

		else if (0 == slice[t])
		{
			f = cx * cy * cz;
			weight[x * size * size + y * size + z] += f;

			f = cx * cy * fz;
			weight[x * size * size + y * size + ((z + 1) % size)] += f;

			f = cx * fy * cz;
			weight[x * size * size + ((y + 1) % size) * size + z] += f;

			f = cx * fy * fz;
			weight[x * size * size + ((y + 1) % size) * size + ((z + 1) % size)] += f;

			f = fx * cy * cz;
			weight[((x + 1) % size) * size * size + y * size + z] += f;

			f = fx * cy * fz;
			weight[((x + 1) % size) * size * size + y * size + ((z + 1) % size)] += f;

			f = fx * fy * cz;
			weight[((x + 1) % size) * size * size + ((y + 1) % size) * size + z] += f;

			f = fx * fy * fz;
			weight[((x + 1) % size) * size * size + ((y + 1) % size) * size + ((z + 1) % size)] += f;
		}
	}
	free(slice);
}

void MakeRotMatrixEuler(double *rot_angle, double rot_matrix[3][3])
{
	// 输入旋转角度rot_angle, 计算出旋转矩阵rot_matrix
	//  rot_angle={alpha, beta, gamma}：以入射X射线方向为z轴负方向，建立右手坐标系，rot_angle表示这个坐标系相对于样品颗粒坐标系的旋转角度，三个角度依次对应于zxz顺规下的欧拉角。
	//  rot_mtrix：由于输入的旋转角是X射线的欧拉角，而rot_matrix代表样品的旋转矩阵，其数值上等于输入欧拉角对应旋转矩阵的逆矩阵。
	double alpha, beta, gamma;
	double s1, c1, s2, c2, s3, c3;

	alpha = rot_angle[0];
	beta = rot_angle[1];
	gamma = rot_angle[2];
	s1 = sin(alpha);
	c1 = cos(alpha);
	s2 = sin(beta);
	c2 = cos(beta);
	s3 = sin(gamma);
	c3 = cos(gamma);

	rot_matrix[0][0] = -s1 * c2 * s3 + c3 * c1;
	rot_matrix[0][1] = -s1 * c2 * c3 - s3 * c1;
	rot_matrix[0][2] = s1 * s2;
	rot_matrix[1][0] = c1 * c2 * s3 + c3 * s1;
	rot_matrix[1][1] = c1 * c2 * c3 - s3 * s1;
	rot_matrix[1][2] = -c1 * s2;
	rot_matrix[2][0] = s2 * s3;
	rot_matrix[2][1] = s2 * c3;
	rot_matrix[2][2] = c2;
}

void PatternCartFloat2Double(float pattern_cart_float[size][size], double pattern_cart_double[size][size])
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			pattern_cart_double[i][j] = pattern_cart_float[i][j];
		}
	}
}

void Volume_1d2Volume_3d(double volume_3d[size][size][size], float *volume_1d, double *weight, int size)
{
	// 根据Merge后的Volume_1d和权重weight计算出三维Volume分布
	memset(volume_3d, 0, sizeof(double) * size * size * size);
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			for (int k = 0; k < size; k++)
			{
				long t = i * size * size + j * size + k;
				if (0 == volume_1d[t])
				{
					volume_3d[i][j][k] = volume_1d[t];
				}
				else
				{
					volume_3d[i][j][k] = volume_1d[t] / weight[t];
				}
			}
		}
	}
}

void Volume3dDouble2FLoat(double volume_3d_double[size][size][size], float volume_3d_float[size][size][size])
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			for (int k = 0; k < size; k++)
			{
				volume_3d_float[i][j][k] = volume_3d_double[i][j][k];
			}
		}
	}
}