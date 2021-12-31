#include <stdlib.h>     /* srand, rand */
#include <math.h>       /* pow */
#include <iostream>
#include <fstream>
#include <string>
#define _USE_MATH_DEFINES
#include <math.h>
#include <random>
using namespace std;

char* load()
// Load a file and return a char pointer where the file was save
{
	ifstream is("lena_256x256.raw", ifstream::binary);
	if (is) {
		is.seekg(0, is.end);
		int length = is.tellg();
		is.seekg(0, is.beg);
		char* buffer = new char[length];
		cout << "Reading " << length << " characters... ";
		is.read(buffer, length);
		//img2 = reinterpret_cast<float*> (buffer);
		is.close();
		return buffer;

	}
	else
	{
		return 0;
	}
}
void save(string filename, char * pointer, float pixel){
// Save an image  file using the pointer where the image is load in memory and the size of the image= width* width * (sizeof pixel))
	ofstream ostrm(filename, ios::binary);
	ostrm.write(pointer, 256 * 256 * (sizeof pixel)); // binary output
	ostrm.close();
}

int main()
{
    const int length=256, width=256;
	int x,y;
    float number;
	float img[length*width], img2[length*width];
	float pixel,pixel2,MSE1,MSE2;
    float sum,sum2=0;
    default_random_engine generator;
    normal_distribution<double> distribution(0.0,0.29); // The varianza that we got from uniform and standar deviation sqrt(varianza)
	for (y = 0; y < width; ++y)
	{
		for (x = 0; x < length; ++x)
		{
            // Uniform distribution
            number=rand() % 100;
            pixel=(number/99)-0.5;
            // Gaussian distribution
            pixel2=distribution(generator);
			img[x + (y*width)] =  pixel;
            img2[x + (y*width)] =  pixel2;
            sum=pow(pixel,2)+sum;
            sum2=pow(pixel2,2)+sum2;
		}
	}
	// Computing and printing (no function) MSE1 and MSE2
    MSE1=sum/(length*width);
    MSE2=sum2/(length*width);
    cout << "MSE= " << MSE1 <<" MSE Gaussian= "<<MSE2<<endl; 
	//Saving noise as files
    char * pointer = reinterpret_cast<char*>(&img);
	save("Uniform.raw",pointer, pixel);
    pointer = reinterpret_cast<char*>(&img2);
    save("Gaussian.raw",pointer, pixel2);
}
