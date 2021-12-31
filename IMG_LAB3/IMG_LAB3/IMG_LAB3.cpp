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

float* matmul(float mat1[256 * 256], float mat2[256 * 256]) {
	float* mul = new float[256 * 256];
	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 256; j++)
		{
			mul[i + 256 * j] = 0;
			for (int k = 0; k < 256; k++)
			{
				mul[i + 256 * j] += mat1[i + 256 * k] * mat2[k + 256 * j];
			}
		}
	}
	return mul;
}

float* transform(float image[256 * 256], float basis[256 * 256]) {

	float* transformed_image = (matmul(basis, image));
	return transformed_image;

}

float* threshold(const float coeff[256 * 256], float tresh) {
	float* thresholded = new float[256 * 256];
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {
			if (abs(coeff[i + 256 * j]) > tresh) {

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

float* transpose(const float matrix[256 * 256]) {
	float* transpose = new float[256 * 256];
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {

			transpose[i + j * 256] = matrix[j + i * 256];

		}
	}
	return transpose;
}

float* transform2d(float image[256 * 256], float basis[256 * 256]) {

	//float* transformed_image = new float[256 * 256];; memory leak

	float* transformed_image = transform(image, basis);
	transformed_image = transpose(transformed_image);
	transformed_image = transform(transformed_image, basis);
	transformed_image = transpose(transformed_image);

	return transformed_image;

}

float* inverse_transform2d(float image[256 * 256], float basis[256 * 256]) {

	//float* transformed_image = new float[256 * 256];; memory leak
	basis = transpose(basis);
	float* transformed_image = transpose(image);
	transformed_image = transform(transformed_image, basis);
	transformed_image = transpose(transformed_image);
	transformed_image = transform(transformed_image, basis);
	return transformed_image;

}

int main()
{
	// 7.1
	const int N =256;
	float pi = 3.14;
	float factor;
	float basis[ N * N];
	

	for (int n = 0; n < N; n++) {
		for (int k = 0; k < N; k++) {
			if (k == 0) {
				factor = sqrt(1.0 / N);
			}
			else {
				factor = sqrt(2.0 / N);
			}
			basis[n + k * N] = factor*cos((pi / N) * (n + 0.5) * k);

		}
	}

	// 7.2
	Store(basis, "basis_vector.raw");
	// 7.3
	float* basis_tr = transpose(basis);
	Store(basis_tr, "basis_vector_tr.raw");
	float* diag = matmul(basis, basis_tr);
	Store(diag, "basis_vector_check.raw");
	// 8.1
	float* lena = load("lena_256x256.raw");
	float* dct_transform = transform(lena, basis);
	Store(dct_transform, "dct_transform.raw");

	float* dct_transform_tr = transpose(dct_transform);
	Store(dct_transform_tr, "dct_transform_tr.raw");
	float* dct_transform_tr_transorm = transform(dct_transform_tr, basis);
	Store(dct_transform_tr_transorm, "dct_transform_tr_transorm.raw");

	//all step are summarize in : 
	dct_transform_tr_transorm = transform2d(lena, basis);
	Store(dct_transform_tr_transorm, "dct_transform_2d.raw");
	dct_transform_tr_transorm = inverse_transform2d(dct_transform_tr_transorm, basis);
	Store(dct_transform_tr_transorm, "dct_inverse_transform_2d.raw");


	//some experiment ..... choose arbitrary treshold 




	//float* transformed = load("transformed_image.raw");
	//threshold(transformed);


	return 0;

}