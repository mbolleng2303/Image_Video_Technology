#include <iostream>
#include <fstream>
#include <string>
#define _USE_MATH_DEFINES
#include <math.h>
using namespace std;

char* load(string filename)
// Load a file and return a char pointer where the file was save
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
void save(string filename, char * pointer, float pixel, int width){
	// Save an image  file using the pointer where the image is load in memory and the size of the image= width* width * (sizeof pixel))
	ofstream ostrm(filename, ios::binary);
	ostrm.write(pointer, width* width * (sizeof pixel)); // binary output
	ostrm.close();
}
void MSEandPSNR(int maxval, float* image1, float* image2, int width){
	//Print the MSE and PSNR comparing two images
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
	const int length=256, width=256;
	int x;
	int y;
	float img[length*width];
	char* charpointer = load("lena_256x256.raw");
	float* floatpointer = reinterpret_cast<float*> (charpointer);
	//Loading pixel values to img raster order;
	float pixel;
	float pixel2;
	for (y = 0; y < width; ++y)	{
		for (x = 0; x < length; ++x){
			pixel =  0.5 + (0.5*cos(x*M_PI / 32)*cos(y*M_PI / 64));
			pixel2=*floatpointer;
			img[x + (y*width)] =  pixel*pixel2;
			floatpointer++;
		}
	}
	floatpointer = reinterpret_cast<float*> (charpointer);
	//Computing and printing MSE and PSNR
	MSEandPSNR(255, floatpointer, (float*) &img, width);
	char * pointer = reinterpret_cast<char*>(&img);
	//Save File Lena*Cos pattern
	save("imageS1.raw",pointer, pixel,width);
}
