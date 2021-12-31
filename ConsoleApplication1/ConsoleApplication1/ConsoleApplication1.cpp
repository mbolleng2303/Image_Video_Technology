
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>

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

void Store(float a[256 * 256], string fname) {
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

	file.write((char*)a, 65536 * 4); //adress in memory

	cout << fname << "writting is good" << endl;
	// closing the file.
	// The reason you need to call close()
	// at the end of the loop is that trying
	// to open a new file without closing the
	// first file will fail.
	file.close();
}

float* matmul(float mat1[8 * 8], float mat2[8 * 8]) {
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
float* transform(float image[8 * 8], float basis[8 * 8]) {

	//float* transformed_image = new float[256 * 256];;

	float* transformed_image = (matmul(basis, image));
	return transformed_image;
	//Store(transformed_image, "transformed_image.raw");
}

float* threshold(const float coeff[8 * 8], float tresh) {
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
float* transpose(const float matrix[8 * 8]) {
	float* transpose = new float[8 * 8];
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++) {

			transpose[i + j * 8] = matrix[j + i * 8];

		}
	}
	return transpose;
}

float* dct_transform2d(float image[8 * 8], float basis[8 * 8]) {

	//float* transformed_image = new float[256 * 256];; memory leak

	float* transformed_image = transform(image, basis);
	transformed_image = transpose(transformed_image);
	transformed_image = transform(transformed_image, basis);
	transformed_image = transpose(transformed_image);

	return transformed_image;
	//Store(transformed_image, "transformed_image.raw");
}
float* inverse_transform2d(float image[8 * 8], float basis[8 * 8]) {

	//float* transformed_image = new float[256 * 256];; memory leak
	basis = transpose(basis);
	float* transformed_image = transpose(image);
	transformed_image = transform(transformed_image, basis);
	transformed_image = transpose(transformed_image);
	transformed_image = transform(transformed_image, basis);
	return transformed_image;
}
float* sliding_window(float image[256 * 256], int l, int k) {
	float* image_as_block = new float[8 * 8];
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++) {
			image_as_block[i + 8 * j] = image[l + i + k + (j * 256)];
		}
	}
	return image_as_block;
}
float* quantize(float dct_transform[8 * 8], int Q[8 * 8]) {
	float* quantized = new float[8 * 8];
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++) {
			quantized[i + 8 * j] = round(dct_transform[i + 8 * j] / Q[i + 8 * j]);
		}
	}
	return quantized;
}
float* inverse_quantize(float img_quantized[8 * 8], float basis[8 * 8], int Q[8 * 8]) {
	float* retrivial = new float[8 * 8];
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++) {
			retrivial[i + 8 * j] = img_quantized[i + 8 * j] * Q[i + 8 * j];
		}
	}
	retrivial = inverse_transform2d(retrivial, basis);
	return retrivial;
}

float* approximate(float image[8 * 8], float basis[8 * 8], int Q[8 * 8]) {

	//float* transformed_image = new float[256 * 256];; memory leak

	float* transformed_image = dct_transform2d(image, basis);
	transformed_image = quantize(transformed_image, Q);

	return transformed_image;
	//Store(transformed_image, "transformed_image.raw");
}


int main()
{
	int Q[8 * 8] = { //given and optimized by hand
		   16,11,10,16,24,40,51,61,
		   12,12,14,19,26,58,60,55,
		   14,13,16,24,40,57,69,56,
		   14,17,22,29,51,87,80,62,
		   18,22,37,56,68,109,103,77,
		   24,35,55,64,81,104,113,92,
		   49,64,78,87,103,121,120,101,
		   72,92,95,98,112,100,103,99 };
	float pi = 3.14;
	float basis[8 * 8];
	for (int n = 0; n < 8; n++) {
		for (int k = 0; k < 8; k++) {

			basis[n + k * 8] = cos((pi / 8) * (n + 0.5) * k);

		}
	}
	Store(basis, "basis_vector.raw");
	float* lena = load("lena_256x256.raw");
	lena = sliding_window(lena, 3, 3);
	Store(lena, "lena_block.raw");
	float* lena_quant = approximate(lena, basis, Q);
	Store(lena_quant, "lena_quantized.raw");
	float* lena_quant_inv = inverse_quantize(lena, basis, Q);
	Store(lena_quant_inv, "lena_quantized_inverse.raw");




}