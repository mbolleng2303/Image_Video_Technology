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
void orthonormal(float* coeff, int length, int width){
    int n,k;
    float (*A)[length] = (float(*)[length]) coeff;
}
int main(){
    const int width=256,length=256; 
    int n,k;
    float coeff[width][length];
    float coeff2[width][length];
    float identity[width][length];
    float factor;
    for (n=0;n<length;++n){
        for (k=0;k<width;++k){
                   
            if(k==0){
                factor= sqrt(1.0/length);
            }else{
                factor= sqrt(2.0/length);
            }
            coeff[n][k]=factor*cos((M_PI /length)*(n+0.5)*k); //generating coefficients
            coeff2[k][n]=coeff[n][k];                         // generating coefficients transpose
        }
    }
    //method to compute the identity matrix
    for (n=0;n<length;++n){
        for (k=0;k<width;++k){
            float sum=0;
            for (int i=0;i<length;++i){
                sum=sum+coeff[n][i]*coeff2[i][k];
            }
            identity[n][k]=sum;
        }
    }
    // Saving identity, DCT coeff and Tranpose
    float element;
    char * pointer = reinterpret_cast<char*>(&coeff);
    save("matrix.raw",pointer,element);
    pointer = reinterpret_cast<char*>(&coeff2);
    save("matrixTransp.raw",pointer,element);
    pointer = reinterpret_cast<char*>(&identity);
    save("identity.raw",pointer,element);
}