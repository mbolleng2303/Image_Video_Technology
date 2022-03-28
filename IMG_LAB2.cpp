#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <chrono>
#include <random>
using namespace std;

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

float MSE(float* img1, float* img2, int size) {
	float mse = 0;
	int x;
	for (x = 0; x < size * size; x++) {

		mse += (img1[x] - img2[x]) * (img1[x] - img2[x]) / (size * size);
	}

	return mse ;
}

float PSNR(int max, float* img1, float* img2 , int size) {

	float mse = MSE(img1, img2, size);

	float psnr = 10 * log10(max * max / mse);

	return psnr;
}

float* gaussian(float mean, float var, int size) {
	// construct a trivial random generator engine from a time-based seed:
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::normal_distribution<float> gaussian_distribution(mean, sqrt(var));
	float* I = new float[size * size];
	std::cout << "random numbers with mean  and variance = " << mean << var << std::endl;
	for (int y = 0; y < size; y++) {
		for (int x = 0; x < size; x++) {
			I[x + y * size] = gaussian_distribution(generator);
		}
	}
	return I;
}
float* uniform(float start, float end, int size) {
	// construct a trivial random generator engine from a time-based seed:
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<float> uniform_distribution(start, end);
	float* I = new float[size * size];
	std::cout << "random numbers with start  and end = " << start << end << std::endl;
	for (int y = 0; y < size; y++) {
		for (int x = 0; x < size; x++) {
			I[x + y * size] = uniform_distribution(generator);
		}
	}
	return I;
}

int main()
{
	int img_size = 256;
	//------------ex4
	//uniform
	float start = -0.5;
	float end = 0.5;
	
	float* U = uniform(start, end, img_size);
	Store(U, "lab2_rand.raw", img_size);

	//gaussian
	float mean = 0;
	float var = 1;
	float* G = gaussian(mean, var, img_size);
	Store(G, "lab2_rand_gaussian.raw", img_size);
	//------------ex5
	
	float* lena = load("lena_256x256.raw");
	float* gaussian_blur = load("lena_256x256_img_blur.raw");
	float mse = MSE(gaussian_blur,lena, img_size);
	std::cout << "mse = " << mse  << std::endl;
	int max = 255;
	float psnr = PSNR(max, lena, gaussian_blur, img_size);
	std::cout << "psnr = " << psnr << std::endl;



	//------------ex6



	int i = 1;
	mse = 0;
	while (mse <=64) {
		float mean = 0;
		float var = i;
		cout << "var" << var << endl;
		float* G = gaussian(mean, var, img_size);

		float* LG = new float[img_size * img_size];
		for (int y = 0; y < img_size; y++) {
			for (int x = 0; x < img_size; x++) {
				LG[x + y * img_size] = lena[x + y * img_size] + G[x + y * img_size];
			}
		}
		Store(LG, "lena_gaussian_65.raw", img_size);
		mse = MSE(LG, lena, img_size);
		std::cout << "mse = " << mse << std::endl;
		max = 255;
		psnr = PSNR(max, lena, LG, img_size);
		std::cout << "psnr = " << psnr << std::endl;
		i += 1;
	}
	//
	float* img_blur = load("lena_gaussian_65_blur.raw");
	mse = MSE(lena, img_blur, img_size);
	std::cout << "mse = " << mse << std::endl;


}
