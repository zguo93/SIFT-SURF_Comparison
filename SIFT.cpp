#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <complex>
#include <vector>
#include <string>
#include <algorithm>
using namespace std;
typedef struct
{
	int x;
	int y;
	int sigma;
	int mainDirec;//main direction of the keypoints
	double offset_x;
	double offset_y;
	double offset_sigma;
	double histogram[36];//36 bins for 360 degrees
	double descriptor[4][4][8];//descriptor with 128 bins. Each 4 x 4 window has 8 bins histogram.
}location;
class Octave
{
public:

	Octave(int row, int col, int number, double** first_img, double blur);
	void blur();
	void diff_Gauss();
	void keypoint_find();
	void keypoint_offset_find(double, double);
	void orientation_assign();
	void descriptor_assembly();


	double** img[6];
	double** blur_img[6];
	double blurfactor;
	double k;//scale for blurfactor
	int index;
	double** gaussian;
	int numberofblur;
	vector <location> keypoints;
	int size_row;
	int size_col;

};
void error(const char*);
int read_pgm_hdr(FILE *, int *, int *);
void **matrix(int, int, int, int, int);
int skipcomment(FILE *fp);
void conv(int row_matrix, int col_matrix, double ** input, double **function, int size_function, double ** output);
void matrix_inverse(double mat[4][4]);
void downsampling(double** img_in, double** img_out, int row, int col);
void output_image(unsigned char**, string, int, int);
void density(Octave first, Octave second, Octave third, Octave fourth, int x_lower, int x_higher, int y_lower, int y_higher);
#define PI 3.14159265



int main(int argc, char **argv)
{
	FILE *inputimage, *outputimage,*inputimage_2;
	int rows, cols, i, j;
	unsigned char **image_in, **image_out,**image_in_2,**image_out_2;
	double **temp_0,**temp_1, **temp_2,**temp_3;// temporary arraries to store the image after downsampling
	double **temp_0_2, **temp_1_2, **temp_2_2, **temp_3_2;
	cout << "SIFT algorithm starts " << endl << " please input the sigma for Gaussian blur kernel (e.g. 0.3) :" << endl;
	double threshold_peak_1, threshold_peak_2, threshold_peak_3, threshold_peak_4;//threshold for peak value in Tylar expansion, each oatave has its only value
	double threshold_determinant_1, threshold_determinant_2, threshold_determinant_3, threshold_determinant_4;
	double blur;//sigma value for gaussian blur kernel
	cin >> blur;
	cout << "threshold for peak value in Tylar expansion for Octave 1(e.g. 0.03 is to remove background points) :" << endl;
	cin >> threshold_peak_1;
	cout << "threshold for peak value in Tylar expansion for Octave 2(e.g. 0) :" << endl;
	cin >> threshold_peak_2;
	cout << "threshold for peak value in Tylar expansion for Octave 3(e.g. 0) :" << endl;
	cin >> threshold_peak_3;
	cout << "threshold for peak value in Tylar expansion for Octave 4(e.g. 0) :" << endl;
	cin >> threshold_peak_4;

	cout << endl;
	cout << "threshold for determinant of Hessian matrix for Octave 1(e.g. 4) :" << endl;
	cin >> threshold_determinant_1;
	cout << "threshold for determinant of Hessian matrix for Octave 2(e.g. 0.02) :" << endl;
	cin >> threshold_determinant_2;
	cout << "threshold for determinant of Hessian matrix for Octave 3(e.g. 0) :" << endl;
	cin >> threshold_determinant_3;
	cout << "threshold for determinant of Hessian matrix for Octave 4(e.g. 0) :" << endl;
	cin >> threshold_determinant_4;



	/* OPEN FILES */
	const char* filename = "IMG_13_1.pgm";
	cout << "the image 1 is " << filename << endl;
	inputimage = fopen(filename, "rb");
	/* READ HEADER */
	if (read_pgm_hdr(inputimage, &rows, &cols) < 0)
		error("not a PGM image or bpp > 8");

	/* OPEN FILES */
	int rows_2, cols_2;
	const char* filename_2 = "IMG_13_4.pgm";
	cout << "the image 2 is " << filename_2 << endl;
	inputimage_2 = fopen(filename_2, "rb");
	/* READ HEADER */
	if (read_pgm_hdr(inputimage_2, &rows_2, &cols_2) < 0)
		error("not a PGM image or bpp > 8");


	/* ALLOCATE ARRAYS */
	image_in = (unsigned char **)matrix(rows , cols, 0, 0, sizeof(char));
	image_out = (unsigned char **)matrix(rows, cols, 0, 0, sizeof(char));
	temp_0 = (double **)matrix(rows, cols, 0, 0, sizeof(double));


	/* READ THE IMAGE */
	for (i = 0; i < rows; i++)
		if (fread(&image_in[i][0], sizeof(char), cols, inputimage) != cols)//load image 1
			error("can't read the image");


	/*initialization*/
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			temp_0[i][j] = image_in[i][j];//convert unsigned char value to floating point number
		}
	}
	clock_t t;
	clock_t total_t;
	total_t = clock();
	t = clock();
	cout << "processing starts " << endl;
	


	Octave first(rows, cols, 1, temp_0, blur);//initialize first octave
	for (int i = 0; i < 5;i++)
		first.blur();
	first.diff_Gauss();
	first.keypoint_find();
	first.keypoint_offset_find(threshold_peak_1, threshold_determinant_1);
	first.orientation_assign();
	first.descriptor_assembly();
	t = clock() - t;
	cout << " Processing first octave has used " << ((float)t) / CLOCKS_PER_SEC << " seconds" << endl;
	cout << "the number of keypoints are " << first.keypoints.size() << endl;
	cout << endl;


	
	if ((rows % 2) == 1)
		rows = (rows - 1) / 2;
	else
		rows = rows / 2;
	if ((cols % 2) == 1)
		cols = (cols - 1) / 2;
	else
		cols = cols / 2;
	temp_1 = (double**)matrix(rows, cols, 0, 0, sizeof(double));
	downsampling(first.img[0],temp_1 , rows, cols);//downsampling
	t = clock();
	Octave second(rows,cols,2,temp_1,blur);//initialize second octave
	for (int i = 0; i < 5; i++)
		second.blur();
	second.diff_Gauss();
	second.keypoint_find();
	second.keypoint_offset_find(threshold_peak_2, threshold_determinant_2);
	second.orientation_assign();
	second.descriptor_assembly();
	t = clock() - t;
	cout << " Processing second octave has used " << ((float)t) / CLOCKS_PER_SEC << " seconds" << endl;
	cout << "the number of keypoints are " << second.keypoints.size() << endl;
	cout << endl;



	if ((rows % 2) == 1)
		rows = (rows - 1) / 2;
	else
		rows = rows / 2;
	if ((cols % 2) == 1)
		cols = (cols - 1) / 2;
	else
		cols = cols / 2;
	temp_2 = (double**)matrix(rows, cols, 0, 0, sizeof(double));
	downsampling(second.img[0], temp_2, rows, cols);
	t = clock();
	Octave third(rows, cols, 3, temp_2, blur);//initialize third octave
	for (int i = 0; i < 5; i++)
		third.blur();
	third.diff_Gauss();
	third.keypoint_find();
	third.keypoint_offset_find(threshold_peak_3, threshold_determinant_3);
	third.orientation_assign();
	third.descriptor_assembly();
	t = clock() - t;
	cout << " Processing third octave has used " << ((float)t) / CLOCKS_PER_SEC << " seconds" << endl;
	cout << "the number of keypoints are " << third.keypoints.size() << endl;
	cout << endl;




	if ((rows % 2) == 1)
		rows = (rows - 1) / 2;
	else
		rows = rows / 2;
	if ((cols % 2) == 1)
		cols = (cols - 1) / 2;
	else
		cols = cols / 2;
	temp_3 = (double**)matrix(rows, cols, 0, 0, sizeof(double));
	downsampling(third.img[0], temp_3, rows, cols);
	t = clock();
	Octave fourth(rows, cols, 4, temp_3, blur);//initialize last octave
	for (int i = 0; i < 5; i++)
		fourth.blur();
	fourth.diff_Gauss();
	fourth.keypoint_find();
	fourth.keypoint_offset_find(threshold_peak_4, threshold_determinant_4);
	fourth.orientation_assign();
	fourth.descriptor_assembly();
	t = clock() - t;
	cout << " Processing fourth octave has used " << ((float)t) / CLOCKS_PER_SEC << " seconds" << endl;
	cout << "the number of keypoints are " << fourth.keypoints.size() << endl;
	cout << endl;

	density(first, second, third, fourth, 50, 100, 120, 220);
	total_t = clock() - total_t;
	cout << "total processing time is " << ((float)total_t) / CLOCKS_PER_SEC << "  seconds" << endl;

	//SIFT for second image
	//
	//

	cout << endl << endl << endl;
	cout << "second image starts to process" << endl;
	total_t = clock();


	image_out_2= (unsigned char **)matrix(rows_2, cols_2, 0, 0, sizeof(char));
	image_in_2 = (unsigned char **)matrix(rows_2, cols_2, 0, 0, sizeof(char));
	temp_0_2 = (double **)matrix(rows_2, cols_2, 0, 0, sizeof(double));

	/* READ THE IMAGE */
	for (i = 0; i < rows_2; i++)
		if (fread(&image_in_2[i][0], sizeof(char), cols_2, inputimage_2) != cols_2)//load image 2
			error("can't read the image");

	/*initialization*/
	for (int i = 0; i < rows_2; i++)
	{
		for (int j = 0; j < cols_2; j++)
		{
			temp_0_2[i][j] = image_in_2[i][j];//convert unsigned char value to floating point number
		}
	}
	t = clock();
	Octave first_2(rows_2, cols_2, 1, temp_0_2, blur);//initialize first octave
	for (int i = 0; i < 5; i++)
		first_2.blur();
	first_2.diff_Gauss();
	first_2.keypoint_find();
	first_2.keypoint_offset_find(threshold_peak_1, threshold_determinant_1);
	first_2.orientation_assign();
	first_2.descriptor_assembly();
	t = clock() - t;
	cout << " Processing fourth octave has used " << ((float)t) / CLOCKS_PER_SEC << " seconds" << endl;
	cout << "the number of keypoints are " << first_2.keypoints.size() << endl;
	cout << endl;

	if ((rows_2 % 2) == 1)
		rows_2 = (rows_2 - 1) / 2;
	else
		rows_2 = rows_2 / 2;
	if ((cols_2 % 2) == 1)
		cols_2 = (cols_2 - 1) / 2;
	else
		cols_2 = cols_2 / 2;
	temp_1_2 = (double**)matrix(rows_2, cols_2, 0, 0, sizeof(double));
	downsampling(first_2.img[0], temp_1_2, rows_2, cols_2);//downsampling
	t = clock();
	Octave second_2(rows_2, cols_2, 2, temp_1_2, blur);//initialize second octave
	for (int i = 0; i < 5; i++)
		second_2.blur();
	second_2.diff_Gauss();
	second_2.keypoint_find();
	second_2.keypoint_offset_find(threshold_peak_2, threshold_determinant_2);
	second_2.orientation_assign();
	second_2.descriptor_assembly();
	t = clock() - t;
	cout << " Processing second octave has used " << ((float)t) / CLOCKS_PER_SEC << " seconds" << endl;
	cout << "the number of keypoints are " << second_2.keypoints.size() << endl;
	cout << endl;

	if ((rows_2 % 2) == 1)
		rows_2 = (rows_2 - 1) / 2;
	else
		rows_2 = rows_2 / 2;
	if ((cols_2 % 2) == 1)
		cols_2 = (cols_2 - 1) / 2;
	else
		cols_2 = cols_2 / 2;
	temp_2_2 = (double**)matrix(rows_2, cols_2, 0, 0, sizeof(double));
	downsampling(second_2.img[0], temp_2_2, rows_2, cols_2);
	t = clock();
	Octave third_2(rows_2, cols_2, 3, temp_2_2, blur);//initialize third octave
	for (int i = 0; i < 5; i++)
		third_2.blur();
	third_2.diff_Gauss();
	third_2.keypoint_find();
	third_2.keypoint_offset_find(threshold_peak_3, threshold_determinant_3);
	third_2.orientation_assign();
	third_2.descriptor_assembly();
	t = clock() - t;
	cout << " Processing third octave has used " << ((float)t) / CLOCKS_PER_SEC << " seconds" << endl;
	cout << "the number of keypoints are " << third_2.keypoints.size() << endl;
	cout << endl;




	if ((rows_2 % 2) == 1)
		rows_2 = (rows_2 - 1) / 2;
	else
		rows_2 = rows_2 / 2;
	if ((cols_2 % 2) == 1)
		cols_2 = (cols_2 - 1) / 2;
	else
		cols_2 = cols_2 / 2;
	temp_3_2 = (double**)matrix(rows_2, cols_2, 0, 0, sizeof(double));
	downsampling(third_2.img[0], temp_3_2, rows_2, cols_2);
	t = clock() - t;
	Octave fourth_2(rows_2, cols_2, 4, temp_3_2,blur);//initialize last octave
	for (int i = 0; i < 5; i++)
		fourth_2.blur();
	fourth_2.diff_Gauss();
	fourth_2.keypoint_find();
	fourth_2.keypoint_offset_find(threshold_peak_4, threshold_determinant_4);
	fourth_2.orientation_assign();
	fourth_2.descriptor_assembly();
	t = clock() - t;
	cout << " Processing fourth octave has used " << ((float)t) / CLOCKS_PER_SEC << " seconds" << endl;
	cout << "the number of keypoints are " << fourth_2.keypoints.size() << endl;
	cout << endl;;

	density(first_2, second_2, third_2, fourth_2, 50, 100, 120, 200);//calculate the density of keypoints
	total_t = clock() - total_t;
	cout << "total processing time is " << ((float)total_t) / CLOCKS_PER_SEC << "  seconds" << endl;



	
	/* WRITE THE IMAGE */
	//the following part is to generate 4 output images. 2 of them are only keypoints locations.
	//the other 2 images combines the image with keypoints. A white dot on images shows a keypoint.
	output_image(image_in, "outputimage.pgm", first.size_row, first.size_col);


	outputimage = fopen("SIFT_points.pgm", "w");
	fprintf(outputimage, "P5\n%d %d\n255\n", first.size_col, first.size_row);
	for (int i = 0; i < first.size_row; i++)
	{
		for (int j = 0; j < first.size_col; j++)
		{
			image_out[i][j] = 0;//convert back to unsigned char

		}
	}
	for (int i = 0; i < first.keypoints.size(); i++)
	{
		image_out[first.keypoints.at(i).x][first.keypoints.at(i).y] = 255;
		//cout <<" first main direction is "<< first.keypoints.at(i).mainDirec << endl;
	}
	for (int i = 0; i < second.keypoints.size(); i++)
	{
		image_out[second.keypoints.at(i).x*2][second.keypoints.at(i).y*2] = 255;
		//cout << " second main direction is " << second.keypoints.at(i).mainDirec << endl;
	}
	for (int i = 0; i < third.keypoints.size(); i++)
	{
		image_out[third.keypoints.at(i).x*4][third.keypoints.at(i).y*4] = 255;
		//cout << " third main direction is " << third.keypoints.at(i).mainDirec << endl;
	}
	for (int i = 0; i < fourth.keypoints.size(); i++)
	{
		image_out[fourth.keypoints.at(i).x*8][fourth.keypoints.at(i).y*8] = 255;
		//cout << " fourth main direction is " << fourth.keypoints.at(i).mainDirec << endl;
	}
	for (i = 0; i < first.size_row; i++)
		if (fwrite(&image_out[i][0], sizeof(char), first.size_col, outputimage) != first.size_col)
			error("can't write the image");
	fclose(outputimage);


	
	outputimage = fopen("SIFT.pgm", "w");
	fprintf(outputimage, "P5\n%d %d\n255\n", first.size_col, first.size_row);
	for (int i = 0; i < first.size_row; i++)
	{
		for (int j = 0; j < first.size_col; j++)
		{
			image_out[i][j] = image_in[i][j];//convert back to unsigned char

		}
	}
	for (int i = 0; i < first.keypoints.size(); i++)
	{
		image_out[first.keypoints.at(i).x][first.keypoints.at(i).y] = 255;
	}
	for (int i = 0; i < second.keypoints.size(); i++)
	{
		image_out[second.keypoints.at(i).x * 2][second.keypoints.at(i).y * 2] = 255;
	}
	for (int i = 0; i < third.keypoints.size(); i++)
	{
		image_out[third.keypoints.at(i).x * 4][third.keypoints.at(i).y * 4] = 255;
	}
	for (int i = 0; i < fourth.keypoints.size(); i++)
	{
		image_out[fourth.keypoints.at(i).x * 8][fourth.keypoints.at(i).y * 8] = 255;
	}
	for (i = 0; i < first.size_row; i++)
		if (fwrite(&image_out[i][0], sizeof(char), first.size_col, outputimage) != first.size_col)
			error("can't write the image");
	fclose(outputimage);



	
	outputimage = fopen("SIFT_points_2.pgm", "w");
	fprintf(outputimage, "P5\n%d %d\n255\n", first_2.size_col, first_2.size_row);
	for (int i = 0; i < first_2.size_row; i++)
	{
		for (int j = 0; j < first_2.size_col; j++)
		{
			image_out_2[i][j] = 0;//convert back to unsigned char

		}
	}
	for (int i = 0; i < first_2.keypoints.size(); i++)
	{
		image_out_2[first_2.keypoints.at(i).x][first_2.keypoints.at(i).y] = 255;
		//cout <<" first main direction is "<< first.keypoints.at(i).mainDirec << endl;
	}
	for (int i = 0; i < second_2.keypoints.size(); i++)
	{
		image_out_2[second_2.keypoints.at(i).x * 2][second_2.keypoints.at(i).y * 2] = 255;
		//cout << " second main direction is " << second.keypoints.at(i).mainDirec << endl;
	}
	for (unsigned i = 0; i < third_2.keypoints.size(); i++)
	{
		image_out_2[third_2.keypoints.at(i).x * 4][third_2.keypoints.at(i).y * 4] = 255;
		//cout << " third main direction is " << third.keypoints.at(i).mainDirec << endl;
	}
	for (unsigned i = 0; i < fourth_2.keypoints.size(); i++)
	{
		image_out_2[fourth_2.keypoints.at(i).x * 8][fourth_2.keypoints.at(i).y * 8] = 255;
		//cout << " fourth main direction is " << fourth.keypoints.at(i).mainDirec << endl;
	}
	for (i = 0; i < first_2.size_row; i++)
		if (fwrite(&image_out_2[i][0], sizeof(char), first_2.size_col, outputimage) != first_2.size_col)
			error("can't write the image");
	fclose(outputimage);


	
	outputimage = fopen("SIFT_2.pgm", "w");
	fprintf(outputimage, "P5\n%d %d\n255\n", first_2.size_col, first_2.size_row);
	for (int i = 0; i < first_2.size_row; i++)
	{
		for (int j = 0; j < first_2.size_col; j++)
		{
			image_out_2[i][j] = image_in_2[i][j];//convert back to unsigned char

		}
	}
	for (unsigned i = 0; i < first_2.keypoints.size(); i++)
	{
		image_out_2[first_2.keypoints.at(i).x][first_2.keypoints.at(i).y] = 255;
	}
	for (unsigned i = 0; i < second_2.keypoints.size(); i++)
	{
		image_out_2[second_2.keypoints.at(i).x * 2][second_2.keypoints.at(i).y * 2] = 255;
	}
	for (unsigned i = 0; i < third_2.keypoints.size(); i++)
	{
		image_out_2[third_2.keypoints.at(i).x * 4][third_2.keypoints.at(i).y * 4] = 255;
	}
	for (unsigned i = 0; i < fourth_2.keypoints.size(); i++)
	{
		image_out_2[fourth_2.keypoints.at(i).x * 8][fourth_2.keypoints.at(i).y * 8] = 255;
	}
	for (i = 0; i < first_2.size_row; i++)
		if (fwrite(&image_out_2[i][0], sizeof(char), first_2.size_col, outputimage) != first_2.size_col)
			error("can't write the image");
	fclose(outputimage);
	
	
	exit(0);
}




void output_image(unsigned char** image_out, string name, int x, int y)
{
	FILE* outputimage;
	outputimage = fopen(name.c_str(), "w");
	fprintf(outputimage, "P5\n%d %d\n255\n", y, x);
	for (int i = 0; i < x; i++)
		if (fwrite(&image_out[i][0], sizeof(char), y, outputimage) !=y)
			error("can't write the image");
	fclose(outputimage);
}

void error(const char *msg)
{
	fprintf(stderr, "%s\n", msg);
	exit(1);
}

int read_pgm_hdr(FILE *fp, int *nrows, int *ncols)
{
	char filetype[3];
	int maxval;

	if (skipcomment(fp) == EOF
		|| fscanf(fp, "%2s\n", filetype) != 1
		|| strcmp(filetype, "P5")
		|| skipcomment(fp) == EOF
		|| fscanf(fp, "%d", ncols) != 1
		|| skipcomment(fp) == EOF
		|| fscanf(fp, "%d", nrows) != 1
		|| skipcomment(fp) == EOF
		|| fscanf(fp, "%d%*c", &maxval) != 1
		|| maxval > 255)
		return(-1);
	else return(0);
}

void **matrix(int nrows, int ncols, int first_row_coord,
	int first_col_coord, int element_size)
{
	void **p;
	int alignment;
	long i;

	if (nrows < 1 || ncols < 1) return(NULL);
	i = nrows*sizeof(void *);
	/* align the addr of the data to be a multiple of sizeof(long double) */
	alignment = i % sizeof(long double);
	if (alignment != 0) alignment = sizeof(long double) - alignment;
	i += nrows*ncols*element_size + alignment;
	if ((p = (void **)malloc((size_t)i)) != NULL)
	{
		/* compute the address of matrix[first_row_coord][0] */
		p[0] = (char *)(p + nrows) + alignment - first_col_coord*element_size;
		for (i = 1; i < nrows; i++)
			/* compute the address of matrix[first_row_coord+i][0] */
			p[i] = (char *)(p[i - 1]) + ncols*element_size;
		/* compute the address of matrix[0][0] */
		p -= first_row_coord;
	}
	return(p);
}


int skipcomment(FILE *fp)
{
	int i;

	if ((i = getc(fp)) == '#')
		while ((i = getc(fp)) != '\n' && i != EOF);
	return(ungetc(i, fp));
}

Octave::Octave(int row, int col, int number, double** first_img,double blur)
{
	//initialize the parameters
	index = number;
	k = sqrt(2);
	if (index == 1)
		blurfactor = blur;
	if (index == 2)
		blurfactor = 2* blur;
	if (index == 3)
		blurfactor = 4* blur;
	if (index == 4)
		blurfactor = 8* blur;
	size_row = row;
	size_col = col;
	numberofblur=0;
	for (int i = 0; i < 6; i++)
	{
		img[i] = (double**)matrix(row, col, 0, 0, sizeof(double));
		blur_img[i] = (double**)matrix(row, col, 0, 0, sizeof(double));
	}
	img[0] = first_img;
}
void Octave::blur()
{
	double sum_gaussian=0;//used for normalization
	if (numberofblur <= 4)
	{
		gaussian = (double**)matrix(7, 7, -3, -3, sizeof(double));
		double blur;
		switch (numberofblur)
		{
		case 0:
			blur = blurfactor; break;
		case 1:
			blur = blurfactor*k; break;
		case 2:
			blur = blurfactor*k*k; break;
		case 3:
			blur = blurfactor*k*k*k; break;
		case 4:
			blur = blurfactor*k*k*k*k; break;
		}
		for (int i = -3; i <= 3; i++)
		{
			for (int j = -3; j <= 3; j++)
			{
				gaussian[i][j] = exp(-(i*i + j*j) / 2 / blur / blur) / (2 * PI*blur*blur);
				//cout << gaussian[i][j] << " ";
				sum_gaussian = sum_gaussian + gaussian[i][j];
			}
			//cout << endl;
		}
		for (int i = -3; i <= 3; i++)
		{
			for (int j = -3; j <= 3; j++)
			{
				gaussian[i][j] = gaussian[i][j] / sum_gaussian;
			}
		}

		std::cout << blur << endl;
		conv(size_row, size_col, img[0], gaussian, 3, img[numberofblur + 1]);
		numberofblur++;
	}
}
void Octave::diff_Gauss()//Get difference of Gaussian images
{
	if (numberofblur != 5)
		std::cout << "blur the image first" << endl;
	else
	{
		for (int i = 1; i <= 4; i++)
		{
			for (int j = 0; j < size_row; j++)
			{
				for (int k = 0; k < size_col; k++)
				{
					blur_img[i][j][k] = img[i][j][k];//store blurred images
					img[i][j][k] = img[i+1][j][k] - img[i][j][k];// store DoG images
				}
			}
		}
	}
}
void Octave::keypoint_find()
{
	int larger,smaller,equal;
	location position;
	for (int i = 2; i < 3; i++)// twice for total 4 imgs
	{
		for (int j = 10; j < size_row-10; j++)//ingore the boundary of img, leaving 10 pixels margin
		{
			for (int k = 10; k < size_col - 10; k++)
			{
				larger = 0;
				smaller = 0;
				equal = 0;
				for (int neighbor_x = -1; neighbor_x <= 1; neighbor_x++)
				{
					for (int neighbor_y = -1; neighbor_y <= 1; neighbor_y++)
					{
						for (int neighbor_z = -1; neighbor_z <= 1; neighbor_z++)
						{
							if (img[i][j][k] == img[i+neighbor_z][j + neighbor_x][k + neighbor_y])
								equal++;
							else if (img[i][j][k] < img[i + neighbor_z][j + neighbor_x][k + neighbor_y])
								smaller++;
							else if (img[i][j][k] > img[i + neighbor_z][j + neighbor_x][k + neighbor_y])
								larger++;
						}
						
					}
				}
				if (smaller==26|| larger== 26)//find the keypoints among 26 neighbors
				{
					position.x = j;//initialize keypoints
					position.y = k;
					position.sigma = i;
					position.offset_sigma = 0;
					position.offset_x = 0;
					position.offset_y = 0;
					position.mainDirec = 0;

					for (int i = 0; i < 36;i++)
						position.histogram[i] =  0 ;
					keypoints.push_back(position);
					//cout << "find " << position.x << "and " << position.y << endl;
				}
					
			}
				
		}
	}
}
void Octave::keypoint_offset_find(double peak, double determinant)
{
	cout << "the initial keypoints " << keypoints.size()<<endl;
	double diff_wrt_x;//difference with respect to row
	double diff_wrt_y;//difference with respect to column
	double diff_wrt_sigma;
	double h[4][4];//3 by 3 matrix of 2nd order derivative
	double peak_value;
	double hessian_m[2][2];
	double trace_hessian;
	double determinant_hessian;
	double temp;//temporary store value
	for (int i = 0; i < keypoints.size(); i++)
	{
		diff_wrt_x = (img[keypoints.at(i).sigma][keypoints.at(i).x + 1][keypoints.at(i).y] - img[keypoints.at(i).sigma][keypoints.at(i).x - 1][keypoints.at(i).y]) / 2;
		diff_wrt_y = (img[keypoints.at(i).sigma][keypoints.at(i).x][keypoints.at(i).y + 1] - img[keypoints.at(i).sigma][keypoints.at(i).x][keypoints.at(i).y - 1]) / 2;
		diff_wrt_sigma = (img[keypoints.at(i).sigma][keypoints.at(i).x][keypoints.at(i).y] - img[keypoints.at(i).sigma-1][keypoints.at(i).x][keypoints.at(i).y]) /(1-sqrt( 2)/2);
		h[2][2] = (img[keypoints.at(i).sigma][keypoints.at(i).x + 1][keypoints.at(i).y] - img[keypoints.at(i).sigma][keypoints.at(i).x][keypoints.at(i).y]) - (img[keypoints.at(i).sigma][keypoints.at(i).x][keypoints.at(i).y] - img[keypoints.at(i).sigma][keypoints.at(i).x - 1][keypoints.at(i).y]);
		h[1][2] = ((img[keypoints.at(i).sigma][keypoints.at(i).x + 1][keypoints.at(i).y + 1] - img[keypoints.at(i).sigma][keypoints.at(i).x + 1][keypoints.at(i).y - 1]) / 2 - (img[keypoints.at(i).sigma][keypoints.at(i).x - 1][keypoints.at(i).y + 1] - img[keypoints.at(i).sigma][keypoints.at(i).x - 1][keypoints.at(i).y - 1]) / 2)/2;
		h[2][3] = ((img[keypoints.at(i).sigma][keypoints.at(i).x + 1][keypoints.at(i).y] - img[keypoints.at(i).sigma-1][keypoints.at(i).x + 1][keypoints.at(i).y]) /(1-sqrt(2)/2)- (img[keypoints.at(i).sigma][keypoints.at(i).x - 1][keypoints.at(i).y] - img[keypoints.at(i).sigma-1][keypoints.at(i).x - 1][keypoints.at(i).y]) / (1-sqrt(2)/2))/2;
		h[2][1] = h[1][2];
		h[1][1] = (img[keypoints.at(i).sigma][keypoints.at(i).x][keypoints.at(i).y + 1] - img[keypoints.at(i).sigma][keypoints.at(i).x][keypoints.at(i).y]) - (img[keypoints.at(i).sigma][keypoints.at(i).x][keypoints.at(i).y] - img[keypoints.at(i).sigma][keypoints.at(i).x][keypoints.at(i).y - 1]);
		h[1][3] = ((img[keypoints.at(i).sigma][keypoints.at(i).x][keypoints.at(i).y + 1] - img[keypoints.at(i).sigma-1][keypoints.at(i).x][keypoints.at(i).y + 1]) / (1 - sqrt(2) / 2) - (img[keypoints.at(i).sigma][keypoints.at(i).x][keypoints.at(i).y - 1] - img[keypoints.at(i).sigma-1][keypoints.at(i).x][keypoints.at(i).y - 1]) / (1 - sqrt(2) / 2)) / 2;
		h[3][2] = h[2][3];
		h[3][1] = h[1][3];
		h[3][3]= ((img[keypoints.at(i).sigma+1][keypoints.at(i).x ][keypoints.at(i).y] - img[keypoints.at(i).sigma][keypoints.at(i).x][keypoints.at(i).y])/(sqrt(2)-1) - (img[keypoints.at(i).sigma][keypoints.at(i).x][keypoints.at(i).y] - img[keypoints.at(i).sigma-1][keypoints.at(i).x][keypoints.at(i).y])/(1-sqrt(2)/2))/(sqrt(2)/2-sqrt(2)/4);
		matrix_inverse(h);
		keypoints.at(i).offset_y = h[1][1] * diff_wrt_y + h[1][2] * diff_wrt_x + h[1][3] * diff_wrt_sigma;
		keypoints.at(i).offset_x = h[2][1] * diff_wrt_y + h[2][2] * diff_wrt_x + h[2][3] * diff_wrt_sigma;
		keypoints.at(i).offset_sigma = h[3][1] * diff_wrt_y + h[3][2] * diff_wrt_x + h[3][3] * diff_wrt_sigma;
		peak_value = img[keypoints.at(i).sigma][keypoints.at(i).x][keypoints.at(i).y] + 0.5*(diff_wrt_x*keypoints.at(i).offset_x + diff_wrt_y*keypoints.at(i).offset_y + diff_wrt_sigma*keypoints.at(i).offset_sigma);
		//cout << "the peak value is " << abs(peak_value) << endl;
		//cout << "offset is " << keypoints.at(i).offset_x <<" " << keypoints.at(i).offset_y << " " << keypoints.at(i).offset_sigma << endl;
		if (keypoints.at(i).offset_x>1 || keypoints.at(i).offset_x < -1 || keypoints.at(i).offset_y>1 || keypoints.at(i).offset_y < -1)
		{
			//cout << "delete " << keypoints.at(i).x << "and " << keypoints.at(i).y << " due to large offset" << endl;
			keypoints.erase(keypoints.begin() + i);
		}
		else if (abs(peak_value)<peak)//delete the keypoints with low contrast locally
		{
			//cout << "delete " << keypoints.at(i).x << "and " << keypoints.at(i).y <<" due to low contrast"<< endl;
			keypoints.erase(keypoints.begin() + i);
			
		}
		else
		{
			hessian_m[0][0] = h[1][1];//Dxx
			hessian_m[0][1] = h[1][2];//Dxy
			hessian_m[1][0] = h[1][2];
			hessian_m[1][1] = h[2][2];//Dyy
			trace_hessian = hessian_m[0][0] + hessian_m[1][1];
			//cout << "trace_hessian " << trace_hessian << endl;
			determinant_hessian = hessian_m[0][0] * hessian_m[1][1] - (hessian_m[0][1] * hessian_m[1][0]);
			//cout << "determinant_hessian " << determinant_hessian << endl;
			if ((trace_hessian*trace_hessian / determinant_hessian > 12.1) || (abs(determinant_hessian)<determinant))//delete the keypoints located on the edge
			{
				//cout << "delete " << keypoints.at(i).x << "and " << keypoints.at(i).y << " located on the edge" << endl;
				keypoints.erase(keypoints.begin() + i);
			}
		}
	}
}

void Octave::orientation_assign()
{
	double** gaussian;
	gaussian = (double**)matrix(51, 51, -25, -25, sizeof(double));
	double blur;
	double** magnitude;
	magnitude = (double**)matrix(51, 51, -25, -25, sizeof(double));
	double** orientation;
	orientation = (double**)matrix(51, 51, -25, -25, sizeof(double));
	int size;
	for (int i = 0; i < keypoints.size(); i++)
	{
		
		blur = 1.5*blurfactor*pow(k,keypoints.at(i).sigma);
		size = blur;
		size = (size - 1) / 2;
		if (keypoints.at(i).x - size - 1>0 && keypoints.at(i).y - size - 1 > 0 && keypoints.at(i).x + size + 1 < size_row &&keypoints.at(i).y + size + 1 < size_col)
		{
			for (int j = -size; j <= size; j++)
			{
				for (int k = -size; k <= size; k++)
				{
					gaussian[j][k] = exp(-(j*j + k*k) / 2 / blur / blur) / (2 * PI*blur*blur);//make gaussian window
					//cout << gaussian[i][j] << " ";
				}
				//cout << endl;
			}
			for (int j = -size; j <= size; j++)//check the orientation in 3x3 region
			{
				for (int k = -size; k <= size; k++)
				{
					magnitude[j][k] = sqrt(pow(blur_img[keypoints.at(i).sigma][keypoints.at(i).x + j + 1][keypoints.at(i).y + k] - blur_img[keypoints.at(i).sigma][keypoints.at(i).x + j - 1][keypoints.at(i).y + k], 2) + pow(blur_img[keypoints.at(i).sigma][keypoints.at(i).x + j][keypoints.at(i).y + k + 1] - blur_img[keypoints.at(i).sigma][keypoints.at(i).x + j][keypoints.at(i).y + k - 1], 2));
					orientation[j][k] = (180 / PI)*atan((blur_img[keypoints.at(i).sigma][keypoints.at(i).x + j + 1][keypoints.at(i).y + k] - blur_img[keypoints.at(i).sigma][keypoints.at(i).x + j - 1][keypoints.at(i).y + k]) / (blur_img[keypoints.at(i).sigma][keypoints.at(i).x + j][keypoints.at(i).y + k + 1] - blur_img[keypoints.at(i).sigma][keypoints.at(i).x + j][keypoints.at(i).y + k - 1]));
					//cout << magnitude[j][k] << "is original and ";
					//make sure the points range from 0 to 360 degrees
					if ((blur_img[keypoints.at(i).sigma][keypoints.at(i).x + j][keypoints.at(i).y + k + 1] - blur_img[keypoints.at(i).sigma][keypoints.at(i).x + j][keypoints.at(i).y + k - 1]) < 0)
						orientation[j][k] = orientation[j][k] + 180;//if the dx is negative, so points are located at 2nd or 3rd quadrant
					if (orientation[j][k] < 0)//for points located at 4th quadrant, change the orientation from [0 -90] to [270 360]
						orientation[j][k] = orientation[j][k] + 360;

					magnitude[j][k] = magnitude[j][k] * gaussian[j][k];
					//cout << magnitude[j][k] << "is final" << endl;
					keypoints.at(i).histogram[(int)(orientation[j][k] / 10)] = keypoints.at(i).histogram[(int)(orientation[j][k] / 10)] + magnitude[j][k];
				}
			}
			double temp = 0;//store the largest value
			for (int num = 0; num < 36; num++)//36 parts, each 10 degrees
			{
				if (keypoints.at(i).histogram[num] >temp)
				{
					temp = keypoints.at(i).histogram[num];//find the largest orientation
					keypoints.at(i).mainDirec = num * 10 + 5;//unit: degree
				}
			}
		}
		else
		{
			keypoints.erase(keypoints.begin() + i);//remove the points too close to boundaries
		}


	}
}
void Octave::descriptor_assembly()
{
	//using 16x16 window to build up sift descriptor
	//Neglect the offset value of keypoints to facilitate the interpolation process. Here the center between 4 pixels will be interpolated as average of 4 neighboring pixels(in correspoding blurred image)
	double** gaussian;
	gaussian = (double**)matrix(18, 18, -9, -9, sizeof(double));
	double** map;
	double** magnitude;
	double** orientation;
	map =(double**)matrix(18, 18, -9, -9, sizeof(double));
	magnitude = (double**)matrix(18, 18, -9, -9, sizeof(double));
	orientation = (double**)matrix(18, 18, -9, -9, sizeof(double));
	for (int j = -8; j <= 7; j++)
	{
		for (int k = -8; k <= 7; k++)
		{
			gaussian[j][k] = exp(-((j+0.5)*(j+0.5) + (k+0.5)*(k+0.5)) / 2 / 8 / 8) / (2 * PI*8*8);//use the sigma equal to one half of the descriptor window, which is 8
			//cout << gaussian[i][j] << " ";
		}
		//cout << endl;
	}
	
	for (int i = 0; i < keypoints.size(); i++)
	{
		for (int j = -9; j <= 8; j++)//interpolate the value and form a matrix
		{
			for (int k = -9; k <= 8; k++)
			{
				map[j][k] = (blur_img[keypoints.at(i).sigma][keypoints.at(i).x+j][keypoints.at(i).y+k]+ blur_img[keypoints.at(i).sigma][keypoints.at(i).x+j+1][keypoints.at(i).y+k]+ blur_img[keypoints.at(i).sigma][keypoints.at(i).x+j+1][keypoints.at(i).y+k+1]+ blur_img[keypoints.at(i).sigma][keypoints.at(i).x+j][keypoints.at(i).y+k+1]) / 4;
			}
		}
		for (int j = -8; j <= 7; j++)//interpolate the value and form a matrix
		{
			for (int k = -8; k <= 7; k++)
			{
				magnitude[j][k] = sqrt(pow((map[j][k + 1] - map[j][k - 1]), 2) + pow((map[j + 1][k] - map[j - 1][k]), 2));

				magnitude[j][k] = magnitude[j][k] * gaussian[j][k];//weight the magnitude based on distance toward center

				orientation[j][k] = (180 / PI)*atan((map[j + 1][k] - map[j - 1][k]) / (map[j][k + 1] - map[j][k - 1]));
				//cout << magnitude[j][k] << "is original and ";
				if ((map[j + 1][k] - map[j - 1][k]) < 0)
					orientation[j][k] = orientation[j][k] + 180;//if the dx is negative, so points are located at 2nd or 3rd quadrant
				if (orientation[j][k] < 0)//for points located at 4th quadrant, change the orientation from [0 -90] to [270 360]
					orientation[j][k] = orientation[j][k] + 360;
				//calculate the relative orientation with respect to main orientation
				orientation[j][k] = orientation[j][k] - keypoints.at(i).mainDirec;
				if (orientation[j][k] < 0)
					orientation[j][k] = orientation[j][k] + 360;
			}
		}
		
		double histogram[8];//8 bin histogram
		for (int x = 0; x < 4; x++)
		{
			for (int y = 0; y < 4; y++)//calculate histogram for each 4x4 window
			{
				for (int number = 0; number < 8; number++)//initialize to zero
					histogram[number] = 0;
				
				for (int j = (x*4) - 8; j < (x*4) - 4; j++)
				{
					for (int k = (y * 4) - 8; k < (y * 4) - 4; k++)
					{
						histogram[(int)(orientation[j][k] / 45)] = histogram[(int)(orientation[j][k] / 45)] + magnitude[j][k];
					}
				}
				
				for (int num = 0; num < 8; num++)
				{
					keypoints.at(i).descriptor[x][y][num] = histogram[num];
				}
			}
		}
		
		//normalize the descriptor
		double sum = 0;
		for (int x = 0; x < 4; x++)
		{
			for (int y = 0; y < 4; y++)
			{
				for (int number = 0; number < 8; number++)
				{
					sum = sum + keypoints.at(i).descriptor[x][y][number];
				}
			}
		}
		for (int x = 0; x <4; x++)
		{
			for (int y = 0; y < 4; y++)
			{
				for (int number = 0; number < 8; number++)
				{
					keypoints.at(i).descriptor[x][y][number] = keypoints.at(i).descriptor[x][y][number] / sum;
					if (keypoints.at(i).descriptor[x][y][number]>0.2)//remove large orientation component
						keypoints.at(i).descriptor[x][y][number] = 0.2;
					
				}

			}
		}
		
		//normalize again
		sum = 0;
		for (int x = 0; x < 4; x++)
		{
			for (int y = 0; y <4; y++)
			{
				for (int number = 0; number < 8; number++)
				{
					sum = sum + keypoints.at(i).descriptor[x][y][number];
				}

			}
		}
		for (int x = 0; x < 4; x++)
		{
			for (int y = 0; y < 4; y++)
			{
				for (int number = 0; number < 8; number++)
				{
					keypoints.at(i).descriptor[x][y][number] = keypoints.at(i).descriptor[x][y][number] / sum;
					//cout << "  " << keypoints.at(i).descriptor[x][y][number] << "  ";
				}

			}
		}
		//cout << endl;
	}

}


void matrix_inverse(double mat[4][4])
{
	//although it requires a 4 by 4 matrix, it only inverse the 3 by 3 part.
	double determinant;
	double temp[4][4];
	determinant = (mat[1][1] * mat[2][2] * mat[3][3]) + (mat[2][1] * mat[3][2] * mat[1][3]) + (mat[3][1] * mat[1][2] * mat[2][3]) - (mat[1][1] * mat[3][2] * mat[2][3]) - (mat[3][1] * mat[2][2] * mat[1][3]) - (mat[2][1] * mat[1][2] * mat[3][3]);
	temp[1][1] = mat[2][2] * mat[3][3] - mat[2][3] * mat[3][2];
	temp[1][2] = mat[1][3] * mat[3][2] - mat[1][2] * mat[3][3];
	temp[1][3] = mat[1][2] * mat[2][3] - mat[1][3] * mat[2][2];
	temp[2][1] = mat[2][3] * mat[3][1] - mat[2][1] * mat[3][3];
	temp[2][2] = mat[1][1] * mat[3][3] - mat[1][3] * mat[3][1];
	temp[2][3] = mat[1][3] * mat[2][1] - mat[1][1] * mat[2][3];
	temp[3][1] = mat[2][1] * mat[3][2] - mat[2][2] * mat[3][1];
	temp[3][2] = mat[1][2] * mat[3][1] - mat[1][1] * mat[3][2];
	temp[3][3] = mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1];
	for (int i = 1; i <= 3; i++)
	{
		for (int j = 1; j <= 3; j++)
			mat[i][j] = temp[i][j]/ determinant;
	}
}
void downsampling(double** img_in, double** img_out, int row, int col)
{
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			img_out[i][j] =img_in[2 * i][2 * j] ;
		}
	}
}


void conv(int row_matrix, int col_matrix, double ** input, double **function, int size_function, double ** output)
{
	// convolution
	double** temp;
	temp = (double**)matrix(row_matrix + 2 * size_function, col_matrix + 2 * size_function, -size_function, -size_function, sizeof(double));
	for (int row = 0; row < row_matrix; row++)
	{
		for (int col = 0; col < col_matrix; col++)
		{
			temp[row][col] = input[row][col];
		}
	}


	for (int row = 0; row < row_matrix; row++)
	{
		for (int col = 0; col < col_matrix; col++)
		{
			double sum = 0;
			for (int i = -size_function; i <= size_function; i++)
			{
				for (int j = -size_function; j <= size_function; j++)
				{

					sum = sum + function[i][j] * temp[row - i][col - j];

				}
			}
			if (sum < 0)
				sum = 0;
			else if (sum>255)
				sum = 255;
			output[row][col] = sum;
		}
	}

}

void density(Octave first, Octave second, Octave third, Octave fourth, int x_lower, int x_higher, int y_lower, int y_higher)
{
	double density = 0;//measure the density of background keypoints
	for (int i = 0; i < first.keypoints.size(); i++)
	{
		if (first.keypoints.at(i).x >= x_lower && first.keypoints.at(i).x < x_higher)
		{
			if (first.keypoints.at(i).y >= y_lower && first.keypoints.at(i).y < y_higher)
				density++;
		}
	}
	for (int i = 0; i < second.keypoints.size(); i++)
	{
		if (second.keypoints.at(i).x >= x_lower && second.keypoints.at(i).x < x_higher)
		{
			if (second.keypoints.at(i).y >= y_lower && second.keypoints.at(i).y < y_higher)
				density++;
		}
	}
	for (int i = 0; i < third.keypoints.size(); i++)
	{
		if (third.keypoints.at(i).x >= x_lower && third.keypoints.at(i).x < x_higher)
		{
			if (third.keypoints.at(i).y >= y_lower && third.keypoints.at(i).y < y_higher)
				density++;
		}
	}
	for (int i = 0; i < fourth.keypoints.size(); i++)
	{
		if (fourth.keypoints.at(i).x >= x_lower && fourth.keypoints.at(i).x < x_higher)
		{
			if (fourth.keypoints.at(i).y >= y_lower && fourth.keypoints.at(i).y < y_higher)
				density++;
		}
	}
	cout << "The number of keypoints is " << density << endl;
	cout << "the density of keypoints at background is " << (density / (x_higher-x_lower) / (y_higher-y_lower)) << endl;
}