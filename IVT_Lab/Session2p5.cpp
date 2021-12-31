#include <stdlib.h>     /* srand, rand */
#include <math.h>       /* pow */
#include <iostream>
#include <fstream>
#include <string>
#define _USE_MATH_DEFINES
#include <math.h>
#include <random>
using namespace std;

char* load(string filename)
{
	ifstream is(filename, ifstream::binary);
	if (is) {
		is.seekg(0, is.end);
		int length = is.tellg();
		is.seekg(0, is.beg);
		char* buffer = new char[length];
		cout << "Reading " << length << " characters... "<<endl;
		is.read(buffer, length);
		is.close();
		return buffer;

	}
	else
	{
		return 0;
	}
}
void save(string filename, char * pointer, float pixel){
	ofstream ostrm(filename, ios::binary);
	ostrm.write(pointer, 256 * 256 * (sizeof pixel)); // binary output
	ostrm.close();
}
void MSEandPSNR(int maxval, float* image1, float* image2, int width){
	int n;
	float err,MSE,PSNR;
	float sum=0;
	for(n=0;n<width*width;n++){
		err=(*image1 - *image2);
		sum=sum+(err*err);
		image1++;
		image2++;
	}
	MSE=sum/(width*width);
	PSNR=10*log10((maxval*maxval)/MSE);
	cout<< "MSE is: " << MSE<<" PSNR is: "<< PSNR <<endl;
}


int main()
{
	char* charpointer = load("lena_256x256.raw");
	float* floatpointer = reinterpret_cast<float*> (charpointer);
    const int length=256, width=256;
	int x,y;
    float number;
	float img[length*width];
	float pixel,pixellena;
    float sum=0,sum2=0;
	float MSE1,MSE2,PSNR1,PSNR2;
	int max=255;
	float std=0.05*max;
    default_random_engine generator;
    normal_distribution<double> distribution(0.0,(std)); // std=5% of the max value
	for (y = 0; y < width; ++y)
	{
		for (x = 0; x < length; ++x)
		{
			pixellena=*floatpointer;
			// Gaussian distribution
			pixel=(distribution(generator));
            img[x + (y*width)] =  pixel+pixellena;
			floatpointer++;
			// 
		}
	}
	//Save Image
    char * pointer;
    pointer = reinterpret_cast<char*>(&img);
    save("Image2lenaNoisy.raw",pointer, pixel);
	floatpointer = reinterpret_cast<float*> (charpointer);
	MSEandPSNR(max, floatpointer, (float*) &img, width);
}