// IMGLAB1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//


#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
using namespace std;

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
float  MSE(float* img1, float* img2, int size) {
	float mse = 0;
	int x;
	for (x = 0; x < size * size; x++)
		mse += (img1[x] - img2[x]) * (img1[x] - img2[x]);

	return mse / (size * size);
}

float PNSR(int max, float* img1, float* img2, int size) {

	float mse = MSE(img1, img2, size);

	float psnr = 10 * log(max * max / mse);

	return psnr;
}
int main()
{
	int img_size = 256;
	//---------------ex1--------
	cout << "Hello" << ' ' << "Wld" << endl;
	//--------------ex2
	float* I = new float[img_size * img_size];
	float pi = 3.14159265359;

	for (int y = 0; y < img_size; y++) {
		for (int x = 0; x < img_size; x++) {
			I[x + y * img_size] = 0.5 + 0.5 * cos(x * pi / 32) * cos(y * pi / 64);
		}
	}
	cout << I << endl;
	cout << "World" << endl;
	Store(I, "lab1.raw", img_size);
	//---------ex3------
	float* op = load("lab1.raw");
	float* lena = load("lena_256x256.raw");
	int x;
	float* lena_new = new float[img_size * img_size];
	for (x = 0; x < img_size * img_size; x++) {
		lena_new[x] = op[x] * lena[x];
	}
	float mse = MSE(lena_new, lena, img_size);
	Store(lena_new, "lena_new.raw", img_size);
	cout << "MSE =" << mse << endl;

}

