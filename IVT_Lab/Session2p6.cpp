#include <stdlib.h>     /* srand, rand */
#include <math.h>       /* pow */
#include <iostream>
#include <fstream>
#include <string>
#define _USE_MATH_DEFINES
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
float kernel(float kernelbox[3][3], float imgbox[3][3]){
	//Method to use a kernel for a 3x3 block
	float sum=0,pixel = 0;
	int k,j;
	for (k = 0; k<3; ++k){
		for (j=0;j<3;++j){
			pixel=pixel +(kernelbox[k][j] * imgbox[k][j]);
		}
	}
	return pixel;
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
	char* charpointer = load("Image2lenaNoisy.raw"); 					//Load the image pointer as char
	float* floatpointer = reinterpret_cast<float*> (charpointer);	// Transform chr to float pointer
    const int length=256, width=256;								// Size of the image 
	float img[length*width], img2[length*width];					// Image memory reservation
	float pixel,pixel2;												// Variables to read pixels
	int x,y;
	//Load pixel values into memory img
    for (y = 0; y < width; ++y)
	{
		for (x = 0; x < length; ++x)
		{
			pixel=*floatpointer;
			img[x + (y*width)] =  pixel;
			floatpointer++;
		}
	}
	//Generating a Gaussian kernel box 
    float Kernelbox[3][3] = {  {0.0751, 0.1238, 0.0751},			// Gaussian 3x3 Kernel sigma=1 
                            {0.1238,0.2042, 0.1238},
                            {0.0751, 0.1238, 0.0751}};
	//Method used to load pixels into blocks and multiply by kernel
	float h_vec1[width]={0};	//vectors to load pixels horizontally
	float h_vec2[width]={0};
	float h_vec3[width]={0};
	for (y = 0; y < width; ++y)
	{
		copy(begin(h_vec2), end(h_vec2), begin(h_vec1));
		if (y<width-1){
			for (x = 0; x < length; ++x){
				h_vec2[x]=img[x + (y*width)];
				h_vec3[x]=img[x + ((y+1)*width)];
			}
		}else{
			float zerovec[width]={0};
			copy(begin(zerovec), end(zerovec), begin(h_vec3));
			for (x = 0; x < length; ++x){
				h_vec2[x]=img[x + (y*width)];
			}
		}
		float imgbox[3][3] = {  {0, 0, 0},
                            	{0, 0, 0},
                            	{0, 0, 0} };
		for (x = 0; x < length; ++x){
			imgbox[0][0]=imgbox[0][1];
			imgbox[1][0]=imgbox[1][1];
			imgbox[2][0]=imgbox[2][1];
			imgbox[0][1]=h_vec1[x];
			imgbox[1][1]=h_vec2[x];
			imgbox[2][1]=h_vec3[x];
			if (x<length-1){
				imgbox[0][2]=h_vec1[x+1];
				imgbox[1][2]=h_vec2[x+1];
				imgbox[2][2]=h_vec2[x+1];
			}else
			{
				imgbox[0][2]=0;
				imgbox[1][2]=0;
				imgbox[2][2]=0;
			}
			img2[x + (y*width)] = kernel(Kernelbox, imgbox); //apply kernel 
		}
	}
	//Save Image
    char * pointer = reinterpret_cast<char*>(&img2);
	save("image2Blur.raw",pointer, pixel);
	//Computing PSNR
	charpointer = load("lena_256x256.raw"); 					//Load the image pointer as char
	float* floatpointer2 = reinterpret_cast<float*> (charpointer);
	cout<<"Before BLUR"<<endl;
	MSEandPSNR(255, floatpointer2, (float*) &img, width);
	cout<<"After BLUR"<<endl;
	floatpointer2 = reinterpret_cast<float*> (charpointer);
	MSEandPSNR(255, floatpointer2, (float*) &img2, width);
}