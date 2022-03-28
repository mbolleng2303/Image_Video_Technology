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

void Store(float* a, string fname, int size) {
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

	file.write((char*)a, size * size * 4); //adress in memory

	cout << fname << "writting is good" << endl;
	// closing the file.
	// The reason you need to call close()
	// at the end of the loop is that trying
	// to open a new file without closing the
	// first file will fail.
	file.close();
}

float* matmul(const float* mat1, const float* mat2, int size) {
	float* mul = new float[size * size];
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			mul[i + size * j] = 0;
			for (int k = 0; k < size; k++)
			{
				mul[i + size * j] += mat1[i + size * k] * mat2[k + size * j];
			}
		}
	}
	return mul;
}

float* transform(const float* image, const float* basis, int size) {

	//float* transformed_image = new float[256 * 256];;

	float* transformed_image = matmul(image, basis, size);
	return transformed_image;
	//Store(transformed_image, "transformed_image.raw");
}

void threshold(float* coeff, const float tresh, int size) {
	//float* thresholded = new float[size * size];
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (abs(coeff[i + size * j]) > tresh) {

				coeff[i + size * j] = coeff[i + size * j];

			}
			else {

				coeff[i + size * j] = 0;
			}
		}
	}
	//return thresholded;

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

float* dct_transform2d(const float* image, const float* basis, int size) {
	float* dct = transform(image, basis, size);
	dct = transpose(dct);
	float* dct2d = transform(dct, basis, size);
	dct2d = transpose(dct2d);
	return dct2d;
}

float* inverse_dct_transform2d(const float* dct2d, const float* basis, int size) {
	const float* basis_t = transpose(basis);
	dct2d = transpose(dct2d);
	float* dct = transform(dct2d, basis_t, size);
	dct = transpose(dct);
	float* image = transform(dct, basis_t, size);
	return image;
}

void quantize(float* dct, const int* Q) {
	for (int i = 0; i < 8; i++)
		for (int j = 0; j < 8; j++)
			dct[i + 8 * j] = round(dct[i + 8 * j] / Q[i + 8 * j]);
}

void inverse_quantize(float* dct, const int* Q) {
	for (int i = 0; i < 8; i++)
		for (int j = 0; j < 8; j++)
			dct[i + 8 * j] = dct[i + 8 * j] * Q[i + 8 * j];
}

float* approximate(const float* image, const float* basis, const int* Q, int size) {
	float* dct = dct_transform2d(image, basis, size);
	threshold(dct, 0, size);
	quantize(dct, Q);
	Store(dct, "quantized.raw", size);
	inverse_quantize(dct, Q);
	return inverse_dct_transform2d(dct, basis, size);
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


float* dicing_window(float* image, int l, int k, int block_size) {
	float* image_as_block = new float[block_size * block_size];
	for (int i = 0; i < block_size; i++) {
		for (int j = 0; j < block_size; j++) {
			image_as_block[i + block_size * j] = image[l + i + ((k + j) * 256)];
		}
	}
	return image_as_block;
}

void insert_block(float* image, const float* block, int k, int l, int block_size) {
	for (int i = 0; i < block_size; i++) {
		for (int j = 0; j < block_size; j++) {

			image[k * block_size + i + (l * block_size + j) * 256] = block[i + block_size * j];
		}
	}
}
void insert_elemen(float* image, float* block, int k, int l, int block_size, int img_size) {
	int img_index = img_size / block_size;
	for (int i = 0; i < block_size; i++) {
		for (int j = 0; j < block_size; j++) {

			image[k  + i * img_index + (l * block_size + j) * 256] = block[i + block_size * j];
		}
	}
}


int main()
{
	const int block_size = 8;
	float pi = 3.14;
	float factor;
	float* basis = new float [block_size * block_size];

	float* Q_mat= new float[8 * 8];
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++) {
			Q_mat[i + 8 * j] = Q[i + 8 * j];
		}
	}

	Store(Q_mat, "Q.raw", 8);
	


	for (int n = 0; n < block_size; n++) {
		for (int k = 0; k < block_size; k++) {
			if (k == 0) {
				factor = sqrt(1.0 / block_size);
			}
			else {
				factor = sqrt(2.0 / block_size);
			}
			basis[n + k * block_size] = factor * cos((pi / block_size) * (n + 0.5) * k);

		}
	}
	// test with 1 block 8x8
	Store(basis, "basis_vector.raw", block_size);

	float* lena = load("lena_256x256.raw");

	float* block = dicing_window(lena, 24, 24, block_size);

	Store(block, "block.raw", block_size);

	float* block_approx = approximate(block, basis, Q, block_size);

	Store(block_approx, "block_approx.raw", block_size);

	//iterative procedure 
	int img_size = 256;
	float* image_retrivial= new float [img_size * img_size];
	
	
	for (int x = 0; x < (img_size / block_size); x++) {
		for (int y = 0; y <( img_size / block_size); y++) {
			
			float* block = dicing_window(lena, x*block_size, y * block_size, block_size);

			float* block_approx = approximate(block, basis, Q, block_size);

			insert_block(image_retrivial, block_approx, x, y, block_size);
			

		}
	}
	Store(image_retrivial, "img_retrivial.raw", img_size);

}