#include <stdlib.h>     /* srand, rand */
#define _USE_MATH_DEFINES /* M_PI*/
#include <math.h>       /* pow */
#include <iostream>
#include <fstream>
using namespace std;
void save(string filename, char * pointer, float pixel){
	ofstream ostrm(filename, ios::binary);
	ostrm.write(pointer, 256 * 256 * (sizeof pixel)); // binary output
	ostrm.close();
}
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
int main(){
    char* charpointer = load("lena_256x256.raw"); 					//Load the image pointer as char
	float* floatpointer = reinterpret_cast<float*> (charpointer);	// Transform chr to float pointer
    const int width=256,length=256; 
    int n,k;
    float coeff[length][width],coeff2[length][width];
    float img[length][width];
    float factor;

    //Loading image in memory and generating DCT coefficients
    for (n=0;n<length;++n){
        for (k=0;k<width;++k){
            if(k==0){
                factor= sqrt(1.0/length);
            }else{
                factor= sqrt(2.0/length);
            }
            coeff[n][k]=factor*cos((M_PI /length)*(n+0.5)*k);
            coeff2[k][n]=coeff[n][k];
            img[n][k]=*floatpointer;
            floatpointer++;
        }
    }
    //Applying DCT 1D
    float imgDCT1[length][width];
    for (n=0;n<length;++n){
        for (k=0;k<width;++k){
            float sum=0.0;
            for (int i=0;i<length;++i){
                sum=sum+img[n][i]*coeff[i][k];
            }
            imgDCT1[n][k]=sum;
        }
    }
    //Applying DCT 2D
    float imgDCT2[length][width];
    for (n=0;n<length;++n){
        for (k=0;k<width;++k){
            float sum=0.0;
            for (int i=0;i<length;++i){
                sum=sum+imgDCT1[i][n]*coeff[i][k];
            }
            imgDCT2[n][k]=sum;
        }
    }
    //Applying Threshold
    float threshold=100;
    for (n=0;n<length;++n){
        for (k=0;k<width;++k){
            if (abs(imgDCT2[n][k])<threshold){
                imgDCT2[n][k]=0;
            }   
        }
    }
    //Saving DCT image 
    float element;
    char * pointer = reinterpret_cast<char*>(&imgDCT2);
    save("lenaDCT2.raw",pointer,element);

    //Inverse DCT
    //In order to reconstruct the image we have to use A inverse =  A transpose--> coeff2
    //Inverse 1D
    float IDCT1[length][width];
    for (n=0;n<length;++n){
        for (k=0;k<width;++k){
            float sum=0.0;
            for (int i=0;i<length;++i){
                sum=sum+imgDCT2[n][i]*coeff2[i][k];
            }
            IDCT1[n][k]=sum;
        }
    }
    //Inverse 2D
    float I[length][width];
    for (n=0;n<length;++n){
        for (k=0;k<width;++k){
            float sum=0.0;
            for (int i=0;i<length;++i){
                sum=sum+IDCT1[i][n]*coeff2[i][k];
            }
            I[n][k]=sum;
        }
    }
    //Saving Image
    pointer = reinterpret_cast<char*>(&I);
    save("IDCT.raw",pointer,element);
}