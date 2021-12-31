#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
//#include <math.h>

using namespace std;

float* load(const char* filename) {
	const int size_in_bytes = 256 * 256 * 4;
	ifstream file(filename, ios::in | ios::binary | ios::ate);// ate:: read at the end 

	if (file.is_open())
	{
		streampos size = file.tellg(); //get the size
		char* memblock = new char[size]; // allocate memory 
		file.seekg(0, ios::beg);// tell the adress of the position 0
		file.read(memblock, size);//
		file.close();

		cout << "the entire file content is in memory with size =" << size << endl;;


		return (float*)memblock;
		//delete[] memblock;
	}
	else cout << "Unable to open file";
	return 0;
}

void Store(float* a, string fname,  int size) {
	// fstream is Stream class to both
	// read and write from/to files.
	// file is object of fstream class
	ofstream file;

	// opening file "fname"
	// in out(write) mode
	// ios::out Open for output operations.
	file.open(fname, ios::binary);

	// If no file is created, then
	// show the error message.
	if (!file)
	{
		cout << "Error in creating file!!!" << endl;

	}

	cout << fname << "File created successfully." << endl;

	file.write((char*)a,size*size * 4); //adress in memory

	cout << fname << "writting is good" << endl;
	// closing the file.
	// The reason you need to call close()
	// at the end of the loop is that trying
	// to open a new file without closing the
	// first file will fail.
	file.close();
}

float* matmul(const float* mat1, const float* mat2) {
	float* mul = new float[8 * 8];
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			mul[i + 8 * j] = 0;
			for (int k = 0; k < 8; k++)
			{
				mul[i + 8 * j] += mat1[i + 8 * k] * mat2[k + 8 * j];
			}
		}
	}
	return mul;
}

float* transform(const float* image, const float* basis) {

	//float* transformed_image = new float[256 * 256];;

	float* transformed_image = matmul(image, basis);
	return transformed_image;
	//Store(transformed_image, "transformed_image.raw");
}

float* threshold(const float* coeff, const float tresh) {
	float* thresholded = new float[8 * 8];
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++) {
			if (abs(coeff[i + 8 * j]) > tresh) {

				thresholded[i] = 0;

			}
			else {

				thresholded[i] = coeff[i];
			}
		}
	}
	return thresholded;

	//Store(img, "threshold_image.raw");
}
float* transpose(const float* matrix) {
	float* transpose = new float[8 * 8];
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++) {

			transpose[i + j * 8] = matrix[j + i * 8];

		}
	}
	return transpose;
}

float* dct_transform2d(const float* image, const float* basis) {
	float* dct = transform(image, basis);
	dct = transpose(dct);
	float* dct2d = transform(dct, basis);
	dct2d = transpose(dct2d);
	return dct2d;
}

float* inverse_dct_transform2d(const float* dct2d, const float* basis) {
	const float* basis_t = transpose(basis);
	dct2d = transpose(dct2d);
	float* dct = transform( dct2d, basis_t);
	dct = transpose(dct);
	float* image = transform(dct, basis_t);
	return image;
}

float* dicing_window(float* image, int l, int k) {
	float* image_as_block = new float[8 * 8];
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++) {
			image_as_block[i + 8 * j] = image[l + i + k + (j * 256)];
		}
	}
	return image_as_block;
}

void quantize(float* dct, const int* Q) {
	for (int i = 0; i < 8; i++)
		for (int j = 0; j < 8; j++)
			dct[i + 8 * j] = dct[i + 8 * j] / Q[i + 8 * j];
}

void inverse_quantize(float* dct, const int* Q) {
	for (int i = 0; i < 8; i++)
		for (int j = 0; j < 8; j++)
			dct[i + 8 * j] = dct[i + 8 * j] * Q[i + 8 * j];
}

float* approximate(const float* image, const float* basis, const int* Q) {
	float* dct = dct_transform2d(image, basis);
	quantize(dct, Q);
	inverse_quantize(dct, Q);
	return inverse_dct_transform2d(dct, basis);
}

int Q[8 * 8] = { //given and optimized by hand
	   16,11,10,16,24,40,51,61,
	   12,12,14,19,26,58,60,55,
	   14,13,16,24,40,57,69,56,
	   14,17,22,29,51,87,80,62,
	   18,22,37,56,68,109,103,77,
	   24,35,55,64,81,104,113,92,
	   49,64,78,87,103,121,120,101,
	   72,92,95,98,112,100,103,99 };

int main()
{
	const int N = 8;
	float pi = 3.14;
	float factor;
	float basis[N * N];


	for (int n = 0; n < N; n++) {
		for (int k = 0; k < N; k++) {
			if (k == 0) {
				factor = sqrt(1.0 / N);
			}
			else {
				factor = sqrt(2.0 / N);
			}
			basis[n + k * N] = factor * cos((pi / N) * (n + 0.5) * k);

		}
	}
	Store(basis, "basis_vector.raw", 8);
	float* lena = load("lena_256x256.raw");

	float* block = dicing_window(lena, 20, 10);
	Store(block, "block.raw", 8);

	float* block_approx = approximate(block, basis, Q);
	Store(block_approx, "block_approx.raw", 8);
	//float*lena_quant_inv = inverse_quantize(block_quant, basis, Q);
	//Store(lena_quant_inv, "lena_quantized_inverse.raw");
}