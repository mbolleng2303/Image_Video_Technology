#include <stdlib.h>     /* srand, rand */
#define _USE_MATH_DEFINES /* M_PI*/
#include <math.h>       /* pow */
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include<sstream> 
#include <map>
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
void Clip(float* img, unsigned char* imgcliped, int width){
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
void Q2RLE(float* QMatrix, int blocksize, float* runlengthvect, int* sizeRunL, int DCVal){
    *runlengthvect=DCVal;
    runlengthvect++;
    if (DCVal<0){
        cout<<"N";
    }
    cout << abs(DCVal);
    float* InterMatrix = QMatrix;
    float Q[8][8];//[blocksize][blocksize];
    int n,k,l,z;
    for (n=0;n<blocksize;n++){
        for (k=0;k<blocksize;k++){
            Q[n][k]=*InterMatrix;
            InterMatrix++;
        }
    }
    z=0;
    l=1;
    n=0;
    k=0;
    int nn=0;
    int kk=1;
    int value;
    // Reading pixels and generating RLE in a diagonal way 
    // there is a problem printing negative values; so we change (-) by (N) means negative symbol
    while (n!=(blocksize-1) || k!=(blocksize-1)){ // run until finish all pixels
        // kk, nn allow to go in diagonal
            if ((k+kk)<0){
                kk=0;
            }else if((n+nn)<0){
                nn=0;
            }
            if((k+kk)>=blocksize){
                k++;
                kk=-1;
                nn=1;
            }else if((n+nn)>=blocksize){
                n++;
                nn=-1;
                kk=1;
            }
            n=n+nn;
            k=k+kk;
            value= Q[n][k];
            if (abs(value)==0){
                if (z==0){
                    cout << z;
                    *runlengthvect=z;
                    runlengthvect++;
                    l++;
                }
                    
                    z++;
            }else{
                if (z==0){
                    if (value<0){
                        cout << "N"; // there is a problem printing negative values; so we change (-) by (N) means negative symbol
                    }
                    cout << abs(value);
                    *runlengthvect=value;
                    runlengthvect++;
                    l++;
                }else{
                    cout << z ;
                    *runlengthvect=z;
                    runlengthvect++;
                    l++;
                    if (value<0){
                        cout << "N"; // there is a problem printing negative values; so we change (-) by (N) means negative symbol
                    }
                    cout << abs(value) ;
                    *runlengthvect=value;
                    runlengthvect++;
                    l++;
                    z=0;
                }
            }
            if (k==7 && n==7 && z!=0){
                cout << z<< endl;
                *runlengthvect=z;
                l++;
                *sizeRunL=*&l;
            }
            if (nn==0){
                kk=-1;
                nn=1;
            } else if(kk==0){
                kk=1;
                nn=-1;
            }
    }

}
void RLE2Q(float* Temporal, float* RLE,int size, int blocksize, int DCval){
    
    int Q[8][8];//[blocksize][blocksize];
    Q[0][0]= DCval+*RLE;
    RLE++;
    int n,k,l;
    int num=0;
    float vector[63];
    l=0;
    for (n=0;n<63;n++){
            if (*RLE!=0){
                vector[n]=*RLE;
                RLE++;
            }else{
                RLE++;
                num=*RLE;
                RLE++;
                for (k=0;k<num;k++){
                     vector[n]=0;
                     n++;
                }
                n--;
            }
    }
    l=0;
    n=0;
    k=0;
    int nn=0;
    int kk=1;
    int value;
    while (n!=(blocksize-1) || k!=(blocksize-1)){
            if ((k+kk)<0){
                kk=0;
            }else if((n+nn)<0){
                nn=0;
            }
            if((k+kk)>=blocksize){
                k++;
                kk=-1;
                nn=1;
            }else if((n+nn)>=blocksize){
                n++;
                nn=-1;
                kk=1;
            }
            n=n+nn;
            k=k+kk;
            Q[n][k]= vector[l] ;
            l++;
            if (nn==0){
                kk=-1;
                nn=1;
            } else if(kk==0){
                kk=1;
                nn=-1;
            }
    }
    for (n=0;n<blocksize;n++){
        for(k=0;k<blocksize;k++){
            *Temporal=Q[n][k];
            Temporal++;
        }
    }
}
void DCT(float* coeff, float* imgblock, float* imgDCT, float* Temporal,int blocksize){
    CosMulti(coeff, imgblock, imgDCT, blocksize);  //Complete 1D DCT
    Transpose(imgDCT,Temporal,blocksize);          // Transpose 1D DCT
    CosMulti(coeff, Temporal, imgDCT, blocksize);  //Complete 2D DCT
}
void encode(float* coeff, float* imgblock, float* imgDCT, float* Temporal, int blocksize,int* sizeRL,float* runlengthvect, int PrevDCVal){
    DCT(coeff, imgblock, imgDCT, Temporal, blocksize);  //Applying DCT 
    Q(imgDCT,Temporal,blocksize);                       //Quantization
    Q2RLE(Temporal,blocksize,runlengthvect, sizeRL, (**&Temporal-PrevDCVal));
}
void decode(float* coeff2,float* imgblock,float* imgDCT, int blocksize, int size, float* RLE, int PrevDCVal){
    float Temporal[8][8];
    RLE2Q(*Temporal, RLE, size,blocksize,PrevDCVal);
    IQ(*Temporal,imgDCT,blocksize);                        //Inverse Quantization
    DCT(coeff2, imgDCT, imgblock, *Temporal, blocksize); //Applying Inverse DCT
}
void approximation(float* coeff, float* coeff2, float* imgblock,float* Temporal, float* imgDCT, int blocksize){
    //encode(coeff, imgblock, imgDCT, Temporal, blocksize);
    
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
    float runlenvect[64][1];
    int sizeRL;
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
    map<const int, int> P;      // P contains all symbols and their frequency 
    map<int, int>::iterator it;
    int val;
    int PrevDCVal=0;
    int maxRL=0;
    for (x=0;  x<width/blocksize; x++){
        for (y=0;  y<length/blocksize; y++){
            ImgBlock(floatpointer,*imgblock, blocksize, x, y, width);
            encode(*coeff, *imgblock, *imgDCT, *Temporal, blocksize, &sizeRL, *runlenvect,PrevDCVal);
            if(maxRL<sizeRL){
                maxRL=sizeRL;
            }
            for (n=0;n<sizeRL;n++){
                val=(int)runlenvect[n][0];
                ++P[val] ;
            }
            decode(*coeff2,*imgblock, *imgDCT, blocksize,sizeRL,*runlenvect,PrevDCVal);
            InsertBlock(*img, *imgblock, blocksize, x, y, width);
            value= *Temporal;
            imgclip[y][x]= *value;
            PrevDCVal=*value;
        }
    }
    float element;
    char * pointer = reinterpret_cast<char*>(&img);
    save("ImgReconstructedRLE.raw",pointer, element, width);
    int totalelements=0;
    int Histsize=P.size();
    for( auto iter = P.begin() ; iter != P.end() ; ++iter )
    {
        cout << "element: " << iter->first << "  frequency: " << iter->second <<endl ;
        totalelements=totalelements+(iter->second);
    }
    float elements[Histsize];
    float Probabilities[Histsize];
    n=0;
    int intvalue;
    float Entropy=0;
    for( auto iter = P.begin() ; iter != P.end() ; ++iter )
    {
        cout << "element: " << iter->first << "  frequency: " << iter->second ;
        elements[n]=(iter->first);
        intvalue=(iter->second);
        Probabilities[n]=(float)((intvalue+0.0)/(totalelements+0.0));
        cout << " P(i) "<<Probabilities[n] <<endl ;
        Entropy=Entropy-(Probabilities[n] * log2(Probabilities[n]));
        n++;
    }
    ///////////////////
    //Printing Histogram uncomment to show symbols as histogram 
    ///////////////////
/*     int plotvalue;
    for (n=0;n<Histsize;n++){
        plotvalue=(int)(Probabilities[n]*700);
        if (elements[n]<0){
            cout<<"N";
        }
        cout<<abs(elements[n]);
        for (k=0;k<plotvalue;k++){
            cout<<"|";
        }
        cout<<endl;
    } */
    // Printing Output values
    cout << "The largest RLE length is: "<< maxRL<< endl;
    cout << "The total number of symbols is: "<< totalelements<< endl;
    cout << "The Number of symbols (Histogram) is: "<< Histsize<< endl;
    cout << "The Average bpsymbol is: "<< Entropy<< endl;
    cout << "The File size is: "<< ((Entropy)*totalelements)/8.0<< " bytes" <<endl;
    // Entropy
}