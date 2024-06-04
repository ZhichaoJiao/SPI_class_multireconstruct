/*20230413*/

#define pi 3.141592653590

extern int size;
extern int n_gamma;
extern int n_r;

void PatternCartFloat2Double(float pattern_cart_float[size][size], double pattern_cart_double[size][size]);

void Volume3dFloat2Double(float volume_3d_float[size][size][size], double volume_3d_double[size][size][size]);

void Volume3dDouble2FLoat(double volume_3d_double[size][size][size], float volume_3d_float[size][size][size]);

void Volume_3d2Volume_1d(double volume_3d[size][size][size], float *volume_1d);

void Volume_1d2Volume_3d(double volume_3d[size][size][size], float *volume_1d, double *weight, int size);

void Normalize(double polar[][n_gamma], int n_r, int n_gamma);

void CorrelationCoefficientPolar_fft_many(const double pattern_polar[][n_gamma], const double reference_polar[][n_gamma], double *correlation_coefficient,
									 int n_r, int n_gamma, int r_min, int r_max);

void Cart2Polar(const double pattern_cart[size][size], double pattern_polar[][n_gamma], double r_min, double r_max, int n_gamma);

void FindBestAngle(const double pattern_polar[][n_gamma], const double reference_polar_matrix[][n_r][n_gamma], double rot_angle_matrix[][3],
				   double best_angle[3], int N_reference, double r_min, double r_max, int n_gamma, double dif_angle_cc[][4], 
				   double *best_angle_cc, double *cc_max);

void MakeAllAngle(double rot_angle_matrix[][3], double step, int *N_angle);

void MakeFineAngle(double best_angle[3], double rot_angle_matrix[][3], double step, double step_fine, int fine_search);

void MakeRotMatrixEuler(double *rot_angle, double rot_matrix[3][3]);

void GenPixels(double *pixels, double lambda, double z_det, double pix_len, int size);

void ReferenceGen(double *rot_angle, double reference_cart[size][size], float *volume_1d, int size, double *pixels);

void PatternMerge(double rot_angle[3], double pattern_cart[size][size], float *volume_1d, double *weight, int size, double *pixels);

void Data2Slice(double pattern_cart[size][size], double *slice, int size);

void Slice2Data(double pattern_cart[size][size], double *slice, int size);

void Hdf5Write3D(hid_t file, int dim[3], char *dataset_name, double ***data);

void Hdf5Write2D(hid_t file, int dim[2], char *dataset_name, double **data);

void Hdf5Write1D(hid_t file, int dim[1], char *dataset_name, double *data);

int min(int x,int y);
