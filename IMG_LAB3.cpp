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

				coeff[i + size * j]= 0;
			}
		}
	}
	//return thresholded;

	//Store(img, "threshold_image.raw");
}

float* transpose(const float* matrix, int size) {
	float* transpose = new float[size * size];
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {

			transpose[i + j * size] = matrix[j + i * size];

		}
	}
	return transpose;
}

float* dct_transform2d(const float* image, const float* basis, int size) {
	float* dct = transform(image, basis, size);
	dct = transpose(dct, size);
	float* dct2d = transform(dct, basis, size);
	dct2d = transpose(dct2d, size);
	return dct2d;
}

float* inverse_dct_transform2d(const float* dct2d, const float* basis, int size) {
	const float* basis_t = transpose(basis, size);
	dct2d = transpose(dct2d, size);
	float* dct = transform(dct2d, basis_t, size);
	dct = transpose(dct, size);
	float* image = transform(dct, basis_t, size);
	return image;
}

float MSE(const float* img1,const float* img2, const int size) {
	float mse = 0;
	int x;
	for (x = 0; x < size * size; x++) {

		mse += (img1[x] - img2[x]) * (img1[x] - img2[x]) / (size * size);
	}

	return mse;
}

float PSNR(const int max, const float* img1, const float* img2, const int size) {

	float mse = MSE(img1, img2, size);

	float psnr = 10 * log10(max * max / mse);

	return psnr;
}

int main()
{
	// 7.1
	const int N =256;
	float pi = 3.14159265359;
	float factor;
	float* basis = new float[ N * N];
	float basis_1D =0;
	
	
	for (int n = 0; n < N; n++) {
		for (int k = 0; k < N; k++) {
			if (k == 0) {
				factor = sqrt(1.0 / N);
			}
			else {
				factor = sqrt( 2.0/ N);
			}
			basis[k + n * N] = factor*cos((pi / N) * (n + 0.5) * k);
			//basis_1D += basis[k + n * N];
		}
	}
	std::cout << "sum = " << basis_1D << std::endl;


	// 7.2
	Store(basis, "basis_vector.raw",N);
	// 7.3
	float* basis_tr = transpose(basis, N);
	Store(basis_tr, "basis_vector_tr.raw",N);
	float* diag = matmul(basis, basis_tr, N);
	Store(diag, "basis_vector_check.raw",N);
	// 8.1
	float* lena = load("lena_256x256.raw");
	float* dct_transform = transform(lena, basis, N);
	Store(dct_transform, "dct_transform.raw",N);

	float* dct_transform_tr = transpose(dct_transform, N);
	Store(dct_transform_tr, "dct_transform_tr.raw",N);

	float* dct_transform_tr_transorm = transform(dct_transform_tr, basis, N);
	Store(dct_transform_tr_transorm, "dct_transform_tr_transorm.raw", N);


	float* dct_transform2 = dct_transform2d(lena, basis, N);
	Store(dct_transform2, "dct_transform_2d.raw",N);
	threshold(dct_transform2, 10 , N);

	Store(dct_transform2, "dct_transform_2d_tresh.raw", N);
	float* inverse_dct_transform2 = inverse_dct_transform2d(dct_transform2, basis, N);
	Store(inverse_dct_transform2, "dct_inverse_transform_2d.raw", N);
	float psnr = PSNR(255, lena, inverse_dct_transform2, N);

	//Store(inverse_dct_transform2, "dct_inverse_transform_2d.raw",N);
	//all step are summarize in : 
	for (float tresh = 0; tresh <50 ; tresh += 5) {
		float* dct_transform2 = dct_transform2d(lena, basis, N);
		//Store(dct_transform2, "dct_transform_2d.raw",N);
		threshold(dct_transform2, tresh, N);

		float* inverse_dct_transform2 = inverse_dct_transform2d(dct_transform2, basis, N);
		float psnr = PSNR(255, lena, inverse_dct_transform2, N);
		//Store(inverse_dct_transform2, "dct_inverse_transform_2d.raw",N);
		std::cout << "tresh" << ' '<< (float)tresh << ' ' <<"psnr" << ' ' << psnr << std::endl;
	}
	return 0;

}