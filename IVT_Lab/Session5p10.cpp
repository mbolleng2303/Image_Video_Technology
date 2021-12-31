#include <stdlib.h>     /* srand, rand */
#define _USE_MATH_DEFINES /* M_PI*/
#include <math.h>       /* pow */
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include<sstream> 
using namespace std;

void save(string filename, char * pointer, float pixel, int imgsize){
	ofstream ostrm(filename, ios::binary);
	ostrm.write(pointer, imgsize * imgsize* (sizeof pixel)); // binary output
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

void Transpose(float* originalMatrix, float* transpMatrix, int imgsize){
    float* img;
    float* transpose=transpMatrix;
    int val;
    int n,k;
    for (n=0;n<imgsize;++n){
        for (k=0;k<imgsize;++k){
            img=originalMatrix + k*imgsize + n;
            *transpose=*img;
            transpose++;
        }
    }
}
void InsertBlock(float* pointer, float* blockpointer, int blocksize, int x, int y, int width){
    float* intpointer=pointer + y*(width*blocksize) + x*blocksize;
    float* intblockpointer=blockpointer;
    int n,k;
    for (n=0;n<blocksize;n++){
        for (k=0;k<blocksize;k++){
            *intpointer=*intblockpointer;
            intblockpointer++;
            intpointer++;
        }
        intpointer=intpointer + width-blocksize;
    }
}
void ImgBlock(float* pointer, float* blockpointer, int blocksize, int x, int y, int width){
    // x and y are the position of the block in the image
    float* intpointer=pointer + y*(width*blocksize) + x*blocksize;
    float* intblockpointer=blockpointer;
    int n,k;
    for (n=0;n<blocksize;n++){
        for (k=0;k<blocksize;k++){
            *intblockpointer=*intpointer;
            intblockpointer++;
            intpointer++;
        }
        intpointer=intpointer + width-blocksize;
    }
}

void CosMulti(float* extcoeff, float* imgblock, float* imgDCTpointer, int blocksize){
    float* intDCTpointer=imgDCTpointer;
    float* img=imgblock;
    float* coeff=extcoeff;
    int n,k;
    for (n=0;n<blocksize;++n){
        for (k=0;k<blocksize;++k){
            float sum=0.0;
            int i;
            float val1,val2;
            for (i=0;i<blocksize;++i){
                val1=*img;
                val2=*coeff;
                sum=sum + val1*val2;
                img++;
                coeff=coeff+blocksize;
            }
            img=img-blocksize;
            coeff=coeff-(i*blocksize)+1;
            *intDCTpointer=sum;
            intDCTpointer++;
        }
        coeff=extcoeff;
        img=img+blocksize;
    }
}
void Q(float* blockimage, float* imgq, int blocksize){
    float* intblockpointer=blockimage;
    float Qat[8][8] = { {16,11,10,16 ,24,40,51,61},	//Q Matrix		
                             {12,12,14,19,26,58,60,55},
                             {14,13,16,24,40,57,69,56},
                             {14,17,22,29,51,87,80,62},
                             {18,22,37,56,68,109,103,77},
                             {24,35,55,64,81,104,113,92},
                             {49,64,78,87,103,121,120,101},
                             {72,92,95,98,112,100,103,99}};
    int n,k;
    for (n=0;n<blocksize;n++){
        for (k=0;k<blocksize;k++){
            *imgq=round(*intblockpointer/Qat[n][k]);
            intblockpointer++;
            *imgq++;
        }
    }
}
void DCT(float* coeff, float* imgblock, float* imgDCT, float* Temporal,int blocksize){
    CosMulti(coeff, imgblock, imgDCT, blocksize);  //Complete 1D DCT
    Transpose(imgDCT,Temporal,blocksize);          // Transpose 1D DCT
    CosMulti(coeff, Temporal, imgDCT, blocksize);  //Complete 2D DCT
}
void encode(float* coeff, float* imgblock, float* imgDCT, float* Temporal, int blocksize){
    DCT(coeff, imgblock, imgDCT, Temporal, blocksize);  //Applying DCT 
    Q(imgDCT,Temporal,blocksize); 
}


int main(){
    
    char* charpointer = load("lena_256x256.raw"); 			//Load the image pointer as char
    float* floatpointer = reinterpret_cast<float*> (charpointer);	// Transform chr to float pointer
    const int width=256, length=256; 
    const int blocksize=8;
    float img[length][width];
    float imgblock[blocksize][blocksize];
    float imgDCT[blocksize][blocksize];
    float Temporal[blocksize][blocksize];
    int n,k;
    float coeff[blocksize][blocksize], coeff2[blocksize][blocksize];
    float factor;
    for (n=0;n<blocksize;++n){
        for (k=0;k<blocksize;++k){
            if(k==0){
                factor= sqrt(1.0/blocksize);
            }else{
                factor= sqrt(2.0/blocksize);
            }
            coeff[n][k]=factor*cos((M_PI /blocksize)*(n+0.5)*k);
            coeff2[k][n]=coeff[n][k];
        }
    }
    const int width2=32,length2=32;
    int x,y;
    float imgclip[length2][width2];
    float* value;
    for (x=0;  x<width/blocksize; x++){
        for (y=0;  y<length/blocksize; y++){
            ImgBlock(floatpointer,*imgblock, blocksize, x, y, width);
            encode(*coeff, *imgblock, *imgDCT, *Temporal, blocksize);
            value= *Temporal;
            imgclip[y][x]= *value;
            
        }
    }
    float element;
    char * pointer = reinterpret_cast<char*>(&imgclip);
    save("ImgDCTerms32x32.raw",pointer, element, width2);
    // Create a text file containing successive differences between quantized DC terms of each block
    int delta[length2*width2];
    vector<int> example;
    ofstream myfile("file.txt");
    int val;
    int l=0;
    for (n=0;n<length2;n++){
        for(k=0;k<width2;k++){
            if (k==0){
                if (n==0){
                    delta[l]=(int) imgclip[n][k];
                    myfile << delta[l] << endl;
                }else{
                    val=(int)(imgclip[n][k]-imgclip[(n-1)][k]);
                    delta[l]=val;
                    myfile << val << endl;
                }
                
            }else{
                val=(int)(imgclip[n][k]-imgclip[n][(k-1)]);
                delta[l]=val;
                myfile << val << endl;
            }
            l++;
        } 
    }
    // Reconstruct Image from differences
    myfile.close();
    ifstream infile;
    infile.open("file.txt");
    int a;
    l=0;
    while (infile >> a)
    {
        delta[l]=a;
        l++;
    }
    float imgQ2[width2][length2];
    l=0;
    for (n=0;n<length2;++n){
        for(k=0;k<width2;k++){
            if (k==0){
                imgQ2[n][k]=delta[l];
                if (n==0){
                    imgQ2[n][k]=delta[l];
                }else{
                    imgQ2[n][k]=(float) delta[l]+imgQ2[(n-1)][k];
                }

            }else{
                imgQ2[n][k]=(float) delta[l]+imgQ2[n][k-1];
            }
            l++;
        }
            
    }
    pointer = reinterpret_cast<char*>(&imgQ2);
    save("ImgRecostructed32X32.raw",pointer,element,width2); 

}