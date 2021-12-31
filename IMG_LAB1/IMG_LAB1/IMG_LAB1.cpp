// IMGLAB1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//



#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
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
	//const int size_in_bytes = 256 * 256 * 4;
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
int main()
{
	//------------ex1
	cout << "Hello" << ' ' << "Wld" << endl;
	//--------------ex2
	float I[256 * 256];
	int x_size = 256;
	int y_size = 256;
	for (int y = 0; y < y_size; y++) {
		for (int x = 0; x < y_size; x++) {
			I[x + y * 256] = 0.5 + 0.5 * cos(x * 3.14 / 32) * cos(y * 3.14 / 64);
		}
	}
	cout << I << endl;
	cout << "World" << endl;
	Store(I, "lab1.raw");
	//---------ex3
	float* op = load("lab1.raw");
	float* lena = load("lena_256x256.raw");
	int x;
	float lena_new[256 * 256];
	//int size = [256 * 256 * 4];
	for (x = 0; x < 256 * 256; x++) {

		lena_new[x] = op[x] * lena[x];


	}
	float mse = MSE(lena_new, lena);
	Store(lena_new, "lena_new.raw");
	cout << "MSE =" << mse << endl;
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