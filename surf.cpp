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
	float offset_x;
	float offset_y;
	float offset_sigma;
	double histogram[36];//36 bins for 360 degrees
	double descriptor[4][4][8];
}location;
typedef struct
{
	int Dxx;
	int Dyy;
	int Dxy;
	double scale;
	double determinant;
}hessian_map;
class Octave
{
public:

	Octave(int row, int col, int number, double** first_img, unsigned char ** original, double threshold);
	void box_blur();
	void find_keypoints();
	void keypoint_offset_find();
	void orientation_assign();
	void descriptor_assembly();


	double** img[5];//0 is original image, 1 to 4 is image for operation
	hessian_map** image[5];//store the hessian matrix result
	double threshold_determinant;
	int index;//for distinguish Octaves
	int size_row;
	int size_col;
	int initial_box_size;// each octave has its box size
	int stepsize_box;
	vector <location> keypoints;


};
void error(const char*);
int read_pgm_hdr(FILE *, int *, int *);
void **matrix(int, int, int, int, int);
int skipcomment(FILE *fp);
void conv(int row_matrix, int col_matrix, double ** input, double **function, int size_function, double ** output);
void Integral(unsigned char**input,double** output, int row, int col);//convert image into an integral image
void matrix_inverse(double mat[4][4]);
void density(Octave first, Octave second, int x_lower, int x_higher, int y_lower, int y_higher);

#define PI 3.14159265



int main(int argc, char **argv)
{
	FILE *inputimage, *outputimage,*inputimage_2,*outputimage_2;
	int rows, cols, i, j,rows_2,cols_2;
	unsigned char **image_in,**image_in_2,**image_out,**image_out_2;
	double **intergal_img,**intergal_img_2;
	double threshold_determinant_1, threshold_determinant_2;
	cout << endl;
	cout << "threshold for determinant of Hessian matrix for Octave 1(e.g. 4*pow(10,3)) :" << endl;
	cin >> threshold_determinant_1;
	cout << "threshold for determinant of Hessian matrix for Octave 2(e.g. 6* pow(10,4)) :" << endl;
	cin >> threshold_determinant_2;

	/* OPEN FILES */
	const char* filename = "IMG_20_1.pgm";
	cout << "the image 1 is " << filename << endl;
	inputimage = fopen(filename, "rb");
	/* READ HEADER */
	if (read_pgm_hdr(inputimage, &rows, &cols) < 0)
		error("not a PGM image or bpp > 8");

	const char* filename_2 = "IMG_20_2.pgm";
	cout << "the image 2 is " << filename_2 << endl;
	inputimage_2 = fopen(filename_2, "rb");
	/* READ HEADER */
	if (read_pgm_hdr(inputimage_2, &rows_2, &cols_2) < 0)
		error("not a PGM image or bpp > 8");

	/* ALLOCATE ARRAYS */
	image_in = (unsigned char **)matrix(rows + 2, cols + 2,-1, -1, sizeof(char));
	intergal_img= (double **)matrix(rows + 2, cols + 2, -1, -1, sizeof(double));
	image_out = (unsigned char **)matrix(rows, cols, 0, 0, sizeof(char));
	image_in_2 = (unsigned char **)matrix(rows_2 + 2, cols_2 + 2, -1, -1, sizeof(char));
	intergal_img_2 = (double **)matrix(rows_2 + 2, cols_2 + 2, -1, -1, sizeof(double));
	image_out_2 = (unsigned char **)matrix(rows_2, cols_2, 0, 0, sizeof(char));

	/* Initialization */
	for (int i = -1; i <= rows; i++)
	{
		for (int j = -1; j <= cols; j++)
		{
			image_in[i][j] = 0;
		}
	}
	for (int i = -1; i <= rows_2; i++)
	{
		for (int j = -1; j <= cols_2; j++)
		{
			image_in_2[i][j] = 0;
		}
	}

	if (image_in == NULL || image_out == NULL) error("can't allocate memory");

	
	/* READ THE IMAGE */
	for (i = 0; i < rows; i++)
		if (fread(&image_in[i][0], sizeof(char), cols, inputimage) != cols)
			error("can't read the image");
	for (i = 0; i < rows_2; i++)
		if (fread(&image_in_2[i][0], sizeof(char), cols_2, inputimage_2) != cols_2)
			error("can't read the image");

	clock_t total_t;
	clock_t t;
	t = clock();
	total_t = clock();
	cout << "processing starts " << endl;

	Integral(image_in, intergal_img, rows, cols);

	Octave first(rows,cols,1, intergal_img, image_in, threshold_determinant_1);
	first.box_blur();
	first.find_keypoints();
	first.keypoint_offset_find();
	first.orientation_assign();
	first.descriptor_assembly();
	t = clock() - t;
	cout << " Processing first octave has used " << ((float)t) / CLOCKS_PER_SEC << " seconds" << endl;
	cout << "first octave: the number of keypoints are " << first.keypoints.size() << endl;
	cout << endl;

	t = clock();
	Octave second(rows, cols, 2, intergal_img, image_in, threshold_determinant_2);
	second.box_blur();
	second.find_keypoints();
	second.keypoint_offset_find();	
	second.orientation_assign();
	second.descriptor_assembly();
	t = clock() - t;
	cout << " Processing second octave has used " << ((float)t) / CLOCKS_PER_SEC << " seconds" << endl;
	cout << "second octave: the number of keypoints are " << second.keypoints.size() << endl;
	cout << endl;


	density(first, second, 50, 100, 120, 200);
	total_t = clock() - total_t;
	cout << "total processing time is " << ((float)total_t) / CLOCKS_PER_SEC << "  seconds" << endl;
	cout << endl << endl << endl;

	cout << "the second image starts to process" << endl;
	total_t = clock();

	Integral(image_in_2, intergal_img_2, rows_2, cols_2);
	t = clock();
	Octave first_2(rows_2, cols_2, 1, intergal_img_2, image_in_2, threshold_determinant_1);
	first_2.box_blur();
	first_2.find_keypoints();
	first_2.keypoint_offset_find();
	first_2.orientation_assign();
	first_2.descriptor_assembly();
	t = clock() - t;
	cout << " Processing first octave has used " << ((float)t) / CLOCKS_PER_SEC << " seconds" << endl;
	cout << "first octave: the number of keypoints are " << first_2.keypoints.size() << endl;
	cout << endl;

	t = clock();
	Octave second_2(rows_2, cols_2, 2, intergal_img_2, image_in_2, threshold_determinant_2);
	second_2.box_blur();
	second_2.find_keypoints();
	second_2.keypoint_offset_find();
	second_2.orientation_assign();
	second_2.descriptor_assembly();
	t = clock() - t;
	cout << " Processing second octave has used " << ((float)t) / CLOCKS_PER_SEC << " seconds" << endl;
	cout << "second octave: the number of keypoints are " << second_2.keypoints.size() << endl;

	density(first_2, second_2, 50, 100, 120, 200);
	total_t = clock() - total_t;
	cout << "total processing time is " << ((float)total_t) / CLOCKS_PER_SEC << "  seconds" << endl;
	cout << endl;


	
	outputimage = fopen("SURF_points.pgm", "w");
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
		//cout << "first octave direction is " << first.keypoints.at(i).mainDirec << endl;

	}
	for (int i = 0; i < second.keypoints.size(); i++)
	{
		image_out[second.keypoints.at(i).x][second.keypoints.at(i).y] = 255;
		//cout <<"second octave direction is "<< second.keypoints.at(i).mainDirec << endl;
	}
	for (i = 0; i < first.size_row; i++)
		if (fwrite(&image_out[i][0], sizeof(char), first.size_col, outputimage) != first.size_col)
			error("can't write the image");
	fclose(outputimage);
	


	
	outputimage = fopen("surf.pgm", "w");
	/* WRITE THE IMAGE */
	fprintf(outputimage, "P5\n%d %d\n255\n", cols, rows);
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
		//cout << "first octave direction is " << first.keypoints.at(i).mainDirec << endl;

	}
	for (int i = 0; i < second.keypoints.size(); i++)
	{
		image_out[second.keypoints.at(i).x][second.keypoints.at(i).y] = 255;
		//cout << "second octave direction is " << second.keypoints.at(i).mainDirec << endl;
	}
	for (i = 0; i < rows; i++)
		if (fwrite(&image_out[i][0], sizeof(char), cols, outputimage) != cols)
			error("can't write the image");


	outputimage_2 = fopen("SURF_points_2.pgm", "w");
	fprintf(outputimage_2, "P5\n%d %d\n255\n", first_2.size_col, first_2.size_row);
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
		//cout << "first octave direction is " << first.keypoints.at(i).mainDirec << endl;

	}
	for (int i = 0; i < second_2.keypoints.size(); i++)
	{
		image_out_2[second_2.keypoints.at(i).x][second_2.keypoints.at(i).y] = 255;
		//cout <<"second octave direction is "<< second.keypoints.at(i).mainDirec << endl;
	}
	for (i = 0; i < first_2.size_row; i++)
		if (fwrite(&image_out_2[i][0], sizeof(char), first_2.size_col, outputimage_2) != first_2.size_col)
			error("can't write the image");
	fclose(outputimage);




	outputimage = fopen("surf_2.pgm", "w");
	/* WRITE THE IMAGE */
	fprintf(outputimage, "P5\n%d %d\n255\n", cols, rows);
	for (int i = 0; i < first_2.size_row; i++)
	{
		for (int j = 0; j < first_2.size_col; j++)
		{
			image_out_2[i][j] = image_in_2[i][j];//convert back to unsigned char

		}
	}
	for (int i = 0; i < first_2.keypoints.size(); i++)
	{
		image_out_2[first_2.keypoints.at(i).x][first_2.keypoints.at(i).y] = 255;
		//cout << "first octave direction is " << first.keypoints.at(i).mainDirec << endl;

	}
	for (int i = 0; i < second_2.keypoints.size(); i++)
	{
		image_out_2[second_2.keypoints.at(i).x][second_2.keypoints.at(i).y] = 255;
		//cout << "second octave direction is " << second.keypoints.at(i).mainDirec << endl;
	}
	for (i = 0; i < rows_2; i++)
		if (fwrite(&image_out_2[i][0], sizeof(char), cols_2, outputimage) != cols_2)
			error("can't write the image");

	
	/* CLOSE FILE & QUIT */
	fclose(outputimage);
	exit(0);
}

void Integral(unsigned char**input, double** output,int row, int col)//compute an integral image for an given input
{
	//create integral image by cumulative addition
	for (int i = -1; i < row+1; i++)
	{
		for (int j = -1; j < col+1; j++)
		{
			output[i][j] = input[i][j];
		}
	}
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			output[i][j] = output[i][j] + output[i][j - 1] + output[i - 1][j] - output[i - 1][j - 1];
		}
	}
}
Octave::Octave(int row, int col, int number, double** first_img, unsigned char ** original,double threshold )
{
	//the constructor of Octave
	//initialize the parameter
	index = number;
	size_row = row;
	size_col = col;
	threshold_determinant = threshold;


	img[0] = (double**)matrix(row, col, 0, 0, sizeof(double));
	image[0] = (hessian_map**)matrix(row, col, 0, 0, sizeof(hessian_map));
	for (int j = 0; j < row; j++)
	{
		for (int k = 0; k < col; k++)
		{
			img[0][j][k] = (double)original[j][k];
			image[0][j][k].Dxx = 0;
			image[0][j][k].Dyy = 0;
			image[0][j][k].Dxy = 0;
			image[0][j][k].determinant = 0;
			image[0][j][k].scale = 0;
		}
	}


	for (int i = 1; i < 5; i++)
	{
		img[i] = (double**)matrix(row, col, 0, 0, sizeof(double));
		image[i]= (hessian_map**)matrix(row, col, 0, 0, sizeof(hessian_map));
		for (int j = 0; j < row; j++)
		{
			for (int k = 0; k < col; k++)
			{
				img[i][j][k] = first_img[j][k];
				image[i][j][k].Dxx = 0;
				image[i][j][k].Dyy = 0;
				image[i][j][k].Dxy = 0;
				image[i][j][k].determinant = 0;
				image[i][j][k].scale = 0;

			}
		}
	}
}


void Octave::box_blur()
{
	if (index == 1)
	{
		initial_box_size = 9; stepsize_box = 6;
	}
	if (index == 2)
	{
		initial_box_size = 15; stepsize_box = 12;
	}
	if (index == 3)
	{
		initial_box_size = 27; stepsize_box = 24;
	}
	int size;//size of box filter
	int length_dxx, height_dxx;
	int length_dyy, height_dyy;
	int length_dxy;
	for (int i = 1; i <= 4; i++)
	{
		size = initial_box_size + (i-1)*(stepsize_box);
		length_dxx = size / 3;
		height_dyy = length_dxx;
		height_dxx = 4 * (size - 9) / 6 + 5;
		length_dyy = height_dxx;
		length_dxy = size / 3;
		for (int j = (size+1)/2; j < (size_row- (size + 1) / 2); j++)
		{
			for (int k = (size + 1) / 2; k < (size_col- (size + 1) / 2); k++)
			{
				image[i][j][k].Dxy = img[i][j - 1][k - 1] + img[i][j - 1 - length_dxy][k - 1 - length_dxy] - (img[i][j - 1 - length_dxy][k - 1] + img[i][j - 1][k - 1 - length_dxy]) + img[i][j][k] + img[i][j + length_dxy][k + length_dxy] - img[i][j][k + length_dxy] - img[i][j + length_dxy][k] + img[i][j - 1][k] + img[i][j - 1 - length_dxy][k + length_dxy] - img[i][j - 1 - length_dxy][k] - img[i][j - 1][k + length_dxy] + img[i][j][k - 1] + img[i][j + length_dxy][k - 1 - length_dxy] - img[i][j][k - 1 - length_dxy] - img[i][j + length_dxy][k - 1 ];
				image[i][j][k].Dxx = -3 * (img[i][j + (length_dxx - 1) / 2][k + (height_dxx - 1) / 2] + img[i][j - (length_dxx - 1) / 2 - 1][k - (height_dxx - 1) / 2 - 1] - img[i][j - (length_dxx - 1) / 2 - 1][k + (height_dxx - 1) / 2] - img[i][j + (length_dxx - 1) / 2][k - (height_dxx - 1) / 2 - 1]) + img[i][j + length_dxx + (length_dxx - 1) / 2][k + (height_dxx - 1) / 2] - img[i][j + length_dxx + (length_dxx - 1) / 2][k + (height_dxx - 1) / 2 - height_dxx] + img[i][j - length_dxx - (length_dxx - 1) / 2 - 1][k - (height_dxx - 1) / 2 - 1] - img[i][j - length_dxx - (length_dxx - 1) / 2 - 1][k - (height_dxx - 1) / 2 - 1+height_dxx];
				image[i][j][k].Dyy = -3 * (img[i][j + (length_dyy - 1) / 2][k + (height_dyy - 1) / 2] + img[i][j - (length_dyy - 1) / 2 - 1][k - (height_dyy - 1) / 2 - 1] - img[i][j - (length_dyy - 1) / 2 - 1][k + (height_dyy - 1) / 2] - img[i][j + (length_dyy - 1) / 2][k - (height_dyy - 1) / 2 - 1]) + img[i][j + (length_dyy - 1) / 2][k + (height_dyy - 1) / 2 + height_dyy] - img[i][j + (length_dyy - 1) / 2 - length_dyy][k + (height_dyy - 1) / 2 + height_dyy] + img[i][j - (length_dyy - 1) / 2 - 1][k - (height_dyy - 1) / 2 - height_dyy - 1] - img[i][j - (length_dyy - 1) / 2 - 1 + length_dyy][k - (height_dyy - 1) / 2 - height_dyy - 1];
				image[i][j][k].determinant = image[i][j][k].Dxx*image[i][j][k].Dyy - 0.9*image[i][j][k].Dxy*0.9*image[i][j][k].Dxy;
				image[i][j][k].scale = 1.2*size / 9;//the same image share the same value
				/*if (i == 4 && index == 1 && j == 498 & k == 500)
				{
					cout << " (height_dxx - 1) / 2 is "<< (height_dxx - 1) / 2<<" img is " << img[i][j + (length_dxx - 1) / 2][k + (height_dxx - 1) / 2] << endl;
					cout << "size length_dxx dyy : "<<size<<"  " <<length_dxx<<"   "<< length_dyy<<endl<<"image[i][j][k].Dxx, Dyy and Dxy are: " << image[i][j][k].Dxx << "  " << image[i][j][k].Dyy << "   " << image[i][j][k].Dxy << endl;

				}*/
					
			}
		}
	}
}
void Octave::find_keypoints()
{
	int larger, smaller, equal;
	location position;
	int size;

	for (int i = 2; i <= 3; i++)
	{
		size = initial_box_size + (i - 1)*(stepsize_box);
		for (int j = (size + 1) / 2+1; j < size_row- (size + 1) / 2 -1 ; j++)
		{
			for (int k = (size + 1) / 2 + 1; k < size_col - (size + 1) / 2 - 1; k++)
			{
				larger = 0;
				smaller = 0;
				equal = 0;
				//if(index==2)
					//cout <<"determinant at "<<j<<" and "<<k<<" is "<< image[i][j][k].determinant << endl;
				for (int neighbor_x = -1; neighbor_x <= 1; neighbor_x++)
				{
					for (int neighbor_y = -1; neighbor_y <= 1; neighbor_y++)
					{
						for (int neighbor_z = -1; neighbor_z <= 1; neighbor_z++)
						{
							if (image[i][j][k].determinant < image[i + neighbor_z][j + neighbor_x][k + neighbor_y].determinant)
								smaller++;
							else if (image[i][j][k].determinant > image[i + neighbor_z][j + neighbor_x][k + neighbor_y].determinant)
								larger++;
						}

					}
				}
				if (smaller == 26 || larger == 26 )
				{
					if (abs(image[i][j][k].determinant) > threshold_determinant)
					{
						//cout << "determinant at " << j << " and " << k << " is " << image[i][j][k].determinant << endl;
						position.x = j;//initialize keypoints
						position.y = k;
						position.sigma = i;
						position.offset_sigma = 0;
						position.offset_x = 0;
						position.offset_y = 0;
						position.mainDirec = 0;
						for (int x = 0; x < 4; x++)
						{
							for (int y = 0; y < 4; y++)
							{
								for (int number = 0; number < 8; number++)
								{
									position.descriptor[x][y][number] = 0;
								}
							}
						}
						for (int x = 0; x < 36; x++)
							position.histogram[x] = 0;
						keypoints.push_back(position);
						//cout << "find " << position.x << " and " << position.y <<" the value is "<<image[i][j][k].determinant << endl;
						//cout <<"the scale is"<< position.sigma << endl;

					}
				}
			}
		}
	}
}
void Octave::keypoint_offset_find()
{
	double diff_wrt_x;
	double diff_wrt_y;
	double diff_wrt_sigma;
	double h[4][4];//3 by 3 matrix of 2nd order derivative
	double peak_value;
	double hessian_m[2][2];
	double trace_hessian;
	double determinant_hessian;

	for (int i = 0; i < keypoints.size(); i++)
	{
		diff_wrt_x = (image[keypoints.at(i).sigma][keypoints.at(i).x + 1][keypoints.at(i).y].determinant - image[keypoints.at(i).sigma][keypoints.at(i).x - 1][keypoints.at(i).y].determinant) / 2;
		diff_wrt_y = (image[keypoints.at(i).sigma][keypoints.at(i).x][keypoints.at(i).y + 1].determinant - image[keypoints.at(i).sigma][keypoints.at(i).x][keypoints.at(i).y - 1].determinant) / 2;
		diff_wrt_sigma = (image[keypoints.at(i).sigma + 1][keypoints.at(i).x][keypoints.at(i).y].determinant - image[keypoints.at(i).sigma - 1][keypoints.at(i).x][keypoints.at(i).y].determinant) / (image[keypoints.at(i).sigma + 1][150][150].scale - image[keypoints.at(i).sigma - 1][150][150].scale);
		h[2][2] = (image[keypoints.at(i).sigma][keypoints.at(i).x + 1][keypoints.at(i).y].determinant - image[keypoints.at(i).sigma][keypoints.at(i).x][keypoints.at(i).y].determinant) - (image[keypoints.at(i).sigma][keypoints.at(i).x][keypoints.at(i).y].determinant - image[keypoints.at(i).sigma][keypoints.at(i).x - 1][keypoints.at(i).y].determinant);
		h[1][2] = ((image[keypoints.at(i).sigma][keypoints.at(i).x + 1][keypoints.at(i).y + 1].determinant - image[keypoints.at(i).sigma][keypoints.at(i).x + 1][keypoints.at(i).y - 1].determinant) / 2 - (image[keypoints.at(i).sigma][keypoints.at(i).x - 1][keypoints.at(i).y + 1].determinant - image[keypoints.at(i).sigma][keypoints.at(i).x - 1][keypoints.at(i).y - 1].determinant) / 2) / 2;
		h[2][3] = ((image[keypoints.at(i).sigma + 1][keypoints.at(i).x + 1][keypoints.at(i).y].determinant - image[keypoints.at(i).sigma - 1][keypoints.at(i).x + 1][keypoints.at(i).y].determinant) / (image[keypoints.at(i).sigma + 1][150][150].scale - image[keypoints.at(i).sigma - 1][150][150].scale) - (image[keypoints.at(i).sigma + 1][keypoints.at(i).x - 1][keypoints.at(i).y].determinant - image[keypoints.at(i).sigma - 1][keypoints.at(i).x - 1][keypoints.at(i).y].determinant) / (image[keypoints.at(i).sigma + 1][150][150].scale - image[keypoints.at(i).sigma - 1][150][150].scale)) / 2;
		h[2][1] = h[1][2];
		h[1][1] = (image[keypoints.at(i).sigma][keypoints.at(i).x][keypoints.at(i).y + 1].determinant - image[keypoints.at(i).sigma][keypoints.at(i).x][keypoints.at(i).y].determinant) - (image[keypoints.at(i).sigma][keypoints.at(i).x][keypoints.at(i).y].determinant - image[keypoints.at(i).sigma][keypoints.at(i).x][keypoints.at(i).y - 1].determinant);
		h[1][3] = ((image[keypoints.at(i).sigma + 1][keypoints.at(i).x][keypoints.at(i).y + 1].determinant - image[keypoints.at(i).sigma - 1][keypoints.at(i).x][keypoints.at(i).y + 1].determinant) / (image[keypoints.at(i).sigma + 1][150][150].scale - image[keypoints.at(i).sigma - 1][150][150].scale) - (image[keypoints.at(i).sigma + 1][keypoints.at(i).x][keypoints.at(i).y - 1].determinant - image[keypoints.at(i).sigma - 1][keypoints.at(i).x][keypoints.at(i).y - 1].determinant) / (image[keypoints.at(i).sigma + 1][150][150].scale - image[keypoints.at(i).sigma - 1][150][150].scale)) / 2;
		h[3][2] = h[2][3];
		h[3][1] = h[1][3];
		h[3][3] = ((image[keypoints.at(i).sigma + 1][keypoints.at(i).x][keypoints.at(i).y].determinant - image[keypoints.at(i).sigma][keypoints.at(i).x][keypoints.at(i).y].determinant) / (image[keypoints.at(i).sigma + 1][150][150].scale - image[keypoints.at(i).sigma][150][150].scale) - (image[keypoints.at(i).sigma][keypoints.at(i).x][keypoints.at(i).y].determinant - image[keypoints.at(i).sigma - 1][keypoints.at(i).x][keypoints.at(i).y].determinant) / (image[keypoints.at(i).sigma][150][150].scale - image[keypoints.at(i).sigma - 1][150][150].scale)) / (0.5*(image[keypoints.at(i).sigma + 1][150][150].scale - image[keypoints.at(i).sigma - 1][150][150].scale));
		matrix_inverse(h);
		keypoints.at(i).offset_x = h[1][1] * diff_wrt_y + h[1][2] * diff_wrt_x + h[1][3] * diff_wrt_sigma;
		keypoints.at(i).offset_y = h[2][1] * diff_wrt_y + h[2][2] * diff_wrt_x + h[2][3] * diff_wrt_sigma;
		keypoints.at(i).offset_sigma = h[3][1] * diff_wrt_y + h[3][2] * diff_wrt_x + h[3][3] * diff_wrt_sigma;
		
		if (keypoints.at(i).offset_x>1 || keypoints.at(i).offset_x < -1 || keypoints.at(i).offset_y>1 || keypoints.at(i).offset_y < -1 || keypoints.at(i).offset_sigma>1 || keypoints.at(i).offset_sigma < -1)
		{
			//cout << "delete " << keypoints.at(i).x << "and " << keypoints.at(i).y << " due to large offset" << endl;
			keypoints.erase(keypoints.begin() + i);
		}
	}
}
void Octave::orientation_assign()
{
	int size;
	int scale;
	double dx;//haar response along y axis
	double dy;//x asix response
	int length;//length for haar wavelet
	double** gaussian;
	gaussian = (double **)matrix(101, 101, -50, -50, sizeof(double));
	double orientation[12] = {0};//divide anlges into 12 parts, each 30 degrees. For reducing computational effect
	double angle;
	double temp;//store for max value in find main orientation
	int temp_orient;
	std::cout <<"for Octave "<<index<< " initial number of keypoints are " << keypoints.size() << endl;
	for (int i = 0; i < keypoints.size(); i++)
	{

		scale = 1.2*((keypoints.at(i).sigma - 1)*stepsize_box + initial_box_size) / 9;
		size = (int)round(6 * scale - 0.5);//get ceiling of floating number
		size = size % 2 == 0 ? size - 1 : size;
		length = (int)round(4 * scale - 0.5);
		length = length % 2 == 0 ? length - 1 : length;//length for calculating the haar transform, rounghly 2*sclae
		dx = 0;//x is vertical axis
		dy = 0;//y is horizontal axis
		temp=0;
		for (int m = 0; m <= 9; m++)
			orientation[m] = 0;
		if ((keypoints.at(i).x - (size - 1) / 2 - (length - 1) / 2>0) && (keypoints.at(i).x + (size - 1) / 2 + (length - 1) / 2 < size_row) && (keypoints.at(i).y - (size - 1) / 2 - (length - 1) / 2 > 0) && (keypoints.at(i).x + (size - 1) / 2 + (length - 1) / 2 < size_col))
		{
			for (int j = -(size - 1) / 2; j <= (size - 1) / 2; j++)
			{
				for (int k = -(size - 1) / 2; k <= (size - 1) / 2; k++)
				{
					if (sqrt(k*k + j*j) <= size)
					{
						gaussian[j][k] = exp(-(j*j + k*k) / 2 / (2 * scale) / (2 * scale)) / (2 * PI*(2 * scale)*(2 * scale));//gaussian function with sigma = 2 scale
						dx = 0;
						dy = 0;
						for (int num = -(length - 1) / 2; num <= (length - 1) / 2; num++)//calculate haar transform on integral image, compute sum
						{
							if (num < 0)
							{
								dx = dx + img[0][keypoints.at(i).x + j + num][keypoints.at(i).y + k];
								dy = dy - img[0][keypoints.at(i).x + j][keypoints.at(i).y + k + num];
							}
							else if (num>0)
							{
								dx = dx - img[0][keypoints.at(i).x + j + num][keypoints.at(i).y + k];
								dy = dy + img[0][keypoints.at(i).x + j][keypoints.at(i).y + k + num];
							}
						}
						//cout << "dx and dy are " << dx << " " << dy << endl;
						if (dy == 0)
							angle = dx > 0 ? 90 : 270;
						else
						{
							angle = (180 / PI)* atan((dx / dy));
						}
						if (dy < 0)
							angle = angle + 180;
						if (angle < 0)
							angle = angle + 360;
						//cout << "position is " << keypoints.at(i).x << " " << keypoints.at(i).y << "      " << j << " " << k  << endl;
						//cout << "angle is " << angle << endl;
						orientation[(int)(angle / 45)] = orientation[(int)(angle / 45)] + gaussian[j][k];

					}

				}
			}
			temp = 0;
			temp_orient = 0;
			for (int x = 0; x < 8; x++)//using wedge of 90 degrees
			{
				//cout << "orientation " << x << " is" << orientation[x]<<"		";//find the orientation with largest weight
				if (x == 0)
				{
					if (temp < (orientation[0] + orientation[7]))
					{
						temp = (orientation[0] + orientation[7]);
						temp_orient = 0;
					}
				}
				else if (temp < (orientation[x - 1] + orientation[x]))
				{
					temp = (orientation[x - 1] + orientation[x]);
					temp_orient = x;
				}
			}
			keypoints.at(i).mainDirec = temp_orient * 45;//Direction is centered at middle of support
			//cout << keypoints.at(i).mainDirec << endl;
		}
		else
		{
			keypoints.erase(keypoints.begin() + i);//erase the points which is too close to boundary
		}

	}
}

void Octave::descriptor_assembly()
{
	int size;//size varies with different octaves and sigmas
	double scale;//scale is the actual sigma value, where sigma stores index of image.
	double sum;//used for normalization
	double dx;//haar response along y axis
	double dy;//x asix response
	double** gaussian;
	double** magnitude;
	double** orientation;
	double** map;
	map = (double**)matrix(261, 261, -130, -130, sizeof(double));//the max possible scale is 13.2, as window is 20scale. so size is around 260.
	gaussian = (double**)matrix(261, 261, -130, -130, sizeof(double));
	magnitude = (double**)matrix(261, 261, -130, -130, sizeof(double));
	orientation = (double**)matrix(261, 261, -130, -130, sizeof(double));


	for (int i = 0; i < keypoints.size(); i++)
	{
		//cout << "sigma is " << keypoints.at(i).sigma << "initial_box size is "<< initial_box_size<<" step is "<<stepsize_box <<endl;
		scale = 1.2*((keypoints.at(i).sigma - 1)*stepsize_box + initial_box_size) / 9;
		size = (int)round(20* scale - 0.5);//get ceiling of floating number
		size = size-1;
		//cout << "size is " << size << endl;
		for (int j = -(size - 1) / 2; j <= (size - 1) / 2; j++)
		{
			for (int k = -(size - 1) / 2; k <= (size - 1) / 2; k++)
			{
				gaussian[j][k] = exp(-(j *j + k*k) / 2 / scale / scale) / (2 * PI * scale * scale);//use the sigma equal to one half of the descriptor window, which is 8

			}
			//cout << endl;
		}
		//cout << "size and x and y are " << (size - 1) / 2 << "  " << keypoints.at(i).x << " " << keypoints.at(i).y << endl;
		if ((keypoints.at(i).x - (size-1)/2 - 1>0 )&& (keypoints.at(i).x + (size - 1) / 2 + 1 < size_row)&&(keypoints.at(i).y - (size - 1) / 2 - 1 > 0) && (keypoints.at(i).y + (size - 1) / 2 + 1 < size_col))
		{
			for (int j = -(size - 1) / 2; j <= (size - 1) / 2; j++)//interpolate the value and form a matrix
			{
				for (int k = -(size - 1) / 2; k <= (size - 1) / 2; k++)
				{
					magnitude[j][k] = sqrt(pow((img[0][keypoints.at(i).x +j][keypoints.at(i).y+k + 1] - img[0][keypoints.at(i).x +j][keypoints.at(i).y+k - 1]), 2) + pow((img[0][keypoints.at(i).x +j + 1][keypoints.at(i).y +k] - img[0][keypoints.at(i).x +j - 1][keypoints.at(i).y+k]), 2));

					magnitude[j][k] = magnitude[j][k] * gaussian[j][k];//weight the magnitude based on distance toward center
					//cout <<"magnitude at "<<j<<" and "<<k <<" is "<<magnitude[j][k] << endl;

					orientation[j][k] = (180 / PI)*atan((img[0][keypoints.at(i).x + j][keypoints.at(i).y + k + 1] - img[0][keypoints.at(i).x + j][keypoints.at(i).y + k - 1]) / (img[0][keypoints.at(i).x + j + 1][keypoints.at(i).y + k] - img[0][keypoints.at(i).x + j - 1][keypoints.at(i).y + k]));
					//cout << magnitude[j][k] << "is original and ";
					if ((img[0][keypoints.at(i).x + j + 1][keypoints.at(i).y + k] - img[0][keypoints.at(i).x + j - 1][keypoints.at(i).y + k]) < 0)
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

					for (int j = (x * (size - 1) / 4 - (size - 1) / 2); j < (x * (size - 1) / 4 - (size - 1) / 4); j++)
					{
						for (int k = (y * (size - 1) / 4 - (size - 1) / 2); k < (y * (size - 1) / 4 - (size - 1) / 4); k++)
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
			sum = 0;// doing normalization
			for (int x = 0; x < 4; x++)
			{
				for (int y = 0; y < 4; y++)//calculate histogram for each 4x4 window
				{
					for (int num = 0; num < 8; num++)
					{
						sum = sum + keypoints.at(i).descriptor[x][y][num];
					}
				}
			}
			for (int x = 0; x < 4; x++)
			{
				for (int y = 0; y < 4; y++)//calculate histogram for each 4x4 window
				{
					for (int num = 0; num < 8; num++)
					{
						keypoints.at(i).descriptor[x][y][num] = keypoints.at(i).descriptor[x][y][num] / sum;
						//cout << " " << keypoints.at(i).descriptor[x][y][num] << " ";
					}
				}
			}
		}else
		{
			//cout << "point was erased" << endl;
			keypoints.erase(keypoints.begin() + i);
		}

	}


}
void matrix_inverse(double mat[4][4])
{
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
			mat[i][j] = temp[i][j] / determinant;
	}
}

void conv(int row_matrix, int col_matrix, double ** input, double **function, int size_function, double ** output)
{
	//convolution
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
void density(Octave first, Octave second,  int x_lower, int x_higher, int y_lower, int y_higher)
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
	cout << "The number of keypoints is " << density << endl;
	cout << "the density of keypoints at background is " << (density / (x_higher - x_lower) / (y_higher - y_lower)) << endl;
}