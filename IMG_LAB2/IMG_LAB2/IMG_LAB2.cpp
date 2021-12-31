// IMGLAB2.cpp : This file contains the 'main' function.Program execution begins and ends there.
//


#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <chrono>
#include <random>
using namespace std;

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

float  MSE(float img1[256 * 256], float img2[256 * 256]) {
	float mse = 0;
	int x;
	for (x = 0; x < 256 * 256; x++)
		mse += (img1[x] - img2[x]) * (img1[x] - img2[x]);

	return mse / (256 * 256);
}

float PNSR(int max, float img1[256 * 256], float img2[256 * 256]) {

	float mse = MSE(img1, img2);

	float psnr = 10 * log(max * max / mse);

	return psnr;
}

float* gaussian(float mean, float var) {
	// construct a trivial random generator engine from a time-based seed:
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::normal_distribution<float> gaussian_distribution(mean, sqrt(var));
	float* I = new float[256 * 256];
	std::cout << "random numbers with mean  and variance = " << mean << var << std::endl;
	for (int y = 0; y < 256; y++) {
		for (int x = 0; x < 256; x++) {
			I[x + y * 256] = gaussian_distribution(generator);
		}
	}
	return (float*)I;
}
float* uniform(float start, float end) {
	// construct a trivial random generator engine from a time-based seed:
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<float> uniform_distribution(start, end);
	float* I = new float[256 * 256];
	std::cout << "random numbers with start  and end = " << start << end << std::endl;
	for (int y = 0; y < 256; y++) {
		for (int x = 0; x < 256; x++) {
			I[x + y * 256] = uniform_distribution(generator);
		}
	}
	return (float*)I;
}

int main()
{
	//------------ex4
	//uniform
	float start = -0.5;
	float end = 0.5;
	float* U = gaussian(start, end);
	Store(U, "lab2_rand.raw");

	//gaussian
	float mean = 0;
	float var = 12;
	float* G = gaussian(mean, var);
	Store(G, "lab2_rand_gaussian.raw");
	//------------ex5
	int x_size = 256;
	int y_size = 256;
	float LG[256 * 256];
	float* lena = load("lena_256x256.raw");
	float* gaussian = load("lab2_rand_gaussian.raw");
	for (int y = 0; y < y_size; y++) {
		for (int x = 0; x < y_size; x++) {
			LG[x + y * 256] = lena[x + y * 256] + gaussian[x + y * 256];
		}
	}
	Store(LG, "lena_gaussian.raw");
	//------------ex6

}



// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file