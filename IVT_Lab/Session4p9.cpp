#include <stdlib.h>     /* srand, rand */
#define _USE_MATH_DEFINES /* M_PI*/
#include <math.h>       /* pow */
#include <iostream>
#include <fstream>
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
    //Method to insert a block of 8x8 into an image 256x256 knowing x,y position
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
    //ImgBlock method to extract a block 8x8 from an image 256x256 knowing x,y
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
    // Method to multiply block of coefficients with an image 
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
void Clip(float* img, unsigned char* imgcliped, int width){
    //Method to clip values form 0-255
    int n,k;
    for (n=0;n<width;n++){
        for(k=0;k<width;k++){
            if (*img>255){
                    *img=255;
            }else if (*img<0){
                    *img=0;
            }
            *imgcliped=(unsigned char) round(*img);
            img++;
            imgcliped++;
        }
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
void IQ(float* blockimage, float* imgq, int blocksize){
    float* intblockpointer=blockimage;
    float Q[8][8] = { {16,11,10,16 ,24,40,51,61},	//Q Matrix		
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
            *imgq=(*intblockpointer) * Q[n][k];
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
    //char * pointer = reinterpret_cast<char*>(imgDCT);
    DCT(coeff, imgblock, imgDCT, Temporal, blocksize);  //Applying DCT 
    //float pixel;
    //save("DCTblock.raw",pointer,pixel,8);
    Q(imgDCT,Temporal,blocksize);                         //Quantization
    //pointer = reinterpret_cast<char*>(Temporal);
    //save("Qlock.raw",pointer,pixel,8);
}

void decode(float* coeff2,float* imgblock,float* Temporal, float* imgDCT, int blocksize){
    IQ(Temporal,imgDCT,blocksize);                        //Inverse Quantization
    //char * pointer = reinterpret_cast<char*>(imgDCT);
    //float pixel;
    //save("IQblock.raw",pointer,pixel,8);
    DCT(coeff2, imgDCT, imgblock, Temporal, blocksize); //Applying Inverse DCT
    //pointer = reinterpret_cast<char*>(imgblock);
    //save("IDCTblock.raw",pointer,pixel,8);
}
void approximation(float* coeff, float* coeff2, float* imgblock,float* Temporal, float* imgDCT, int blocksize){
    encode(coeff, imgblock, imgDCT, Temporal, blocksize);
    decode(coeff2,imgblock,Temporal, imgDCT, blocksize);
}




int main(){
    char* charpointer = load("lena_256x256.raw"); 					//Load the image pointer as char
	float* floatpointer = reinterpret_cast<float*> (charpointer);	// Transform chr to float pointer
    const int width=256,length=256; 
    const int blocksize=8;
    float img[length][width];
    unsigned char imgclip[length][width];
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
    // Procedure to code and decode image block by block x,y are the position of the block
    int x,y;
    for (x=0; x<width/blocksize;x++){
        for (y=0; y<length/blocksize;y++){
            
            
            ImgBlock(floatpointer,*imgblock,blocksize, x, y, width);
            //float pixel;
            //char * pointer = reinterpret_cast<char*>(&imgblock);
            //save("imgblock.raw",pointer,pixel,8);
            approximation(*coeff, *coeff2, *imgblock,*Temporal, *imgDCT, blocksize);
            InsertBlock(*img, *imgblock, blocksize, x, y, width);
        }
    }
    float element;
    char * pointer = reinterpret_cast<char*>(&img);
    save("lenajpeg.raw",pointer, element, width);
    Clip(*img, *imgclip,width);
    unsigned char element2;
    pointer = reinterpret_cast<char*>(&imgclip);
    save("lena8bpp.raw",pointer, element2, width);
}