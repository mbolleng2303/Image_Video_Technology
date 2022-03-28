#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
//#include <math.h>
#include <string>
#include <map>
#include <iterator>
#include<sstream> 

using namespace std;

float* load(const char* filename) {
	const int size_in_bytes = 256 * 256 * 4;
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

void Store(float* a, string fname, int size) {
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

	file.write((char*)a, size * size * 4); //adress in memory

	cout << fname << "writting is good" << endl;
	// closing the file.
	// The reason you need to call close()
	// at the end of the loop is that trying
	// to open a new file without closing the
	// first file will fail.
	file.close();
}
float MSE(const float* img1, const float* img2, const int size) {
	float mse = 0;
	int x;
	for (x = 0; x < size * size; x++) {

		mse += (img1[x] - img2[x]) * (img1[x] - img2[x]) / (size * size);
	}

	return mse;
}
float PSNR(const int max, const float* img1, const float* img2, const int size) {

	float mse = MSE(img1, img2, size);

	float psnr = 10 * log10(max * max / mse);

	return psnr;
}
float* matmul(const float* mat1, const float* mat2, int size) {
	float* mul = new float[size * size];
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			mul[i + size * j] = 0;
			for (int k = 0; k < size; k++)
			{
				mul[i + size * j] += mat1[i + size * k] * mat2[k + size * j];
			}
		}
	}
	return mul;
}
float* transform(const float* image, const float* basis, int size) {

	//float* transformed_image = new float[256 * 256];;

	float* transformed_image = matmul(image, basis, size);
	return transformed_image;
	//Store(transformed_image, "transformed_image.raw");
}
void threshold(float* coeff, const float tresh, int size) {
	//float* thresholded = new float[size * size];
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (abs(coeff[i + size * j]) > tresh) {

				coeff[i + size * j] = coeff[i + size * j];

			}
			else {

				coeff[i + size * j] = 0;
			}
		}
	}
	//return thresholded;

	//Store(img, "threshold_image.raw");
}
float* transpose(const float* matrix, int size) {
	float* transpose = new float[size * size];
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {

			transpose[i + j * size] = matrix[j + i * size];

		}
	}
	return transpose;
}

float* dct_transform2d(const float* image, const float* basis, int size) {
	float* dct = transform(image, basis, size);
	dct = transpose(dct, size);
	float* dct2d = transform(dct, basis, size);
	dct2d = transpose(dct2d, size);
	return dct2d;
}

float* inverse_dct_transform2d(const float* dct2d, const float* basis, int size) {
	const float* basis_t = transpose(basis, size);
	dct2d = transpose(dct2d, size);
	float* dct = transform(dct2d, basis_t, size);
	dct = transpose(dct, size);
	float* image = transform(dct, basis_t, size);
	return image;
}

void quantize(float* dct, const int* Q) {
	for (int i = 0; i < 8; i++)
		for (int j = 0; j < 8; j++)
			dct[i + 8 * j] = round(dct[i + 8 * j] / Q[i + 8 * j]);
}

void inverse_quantize(float* dct, const int* Q) {
	for (int i = 0; i < 8; i++)
		for (int j = 0; j < 8; j++)
			dct[i + 8 * j] = dct[i + 8 * j] * Q[i + 8 * j];
}

float* encode(const float* image, const float* basis, const int* Q, int size) {
	float* dct = dct_transform2d(image, basis, size);
	threshold(dct, 10, size);
	quantize(dct, Q);
	return dct;
}

float* decode(float* dct, const float* basis, const int* Q, int size) {
	inverse_quantize(dct, Q);
	return inverse_dct_transform2d(dct, basis, size);
}

float* approximate(const float* image, const float* basis, const int* Q, int size) {

	float* dct = encode(image, basis, Q, size);
	float* img_approx = decode(dct, basis, Q, size);
	return img_approx;
}
const float* get_basis(int size) {
	const int block_size = size;
	float pi = 3.14;
	float factor;
	float* basis = new float[block_size * block_size];
	for (int n = 0; n < block_size; n++) {
		for (int k = 0; k < block_size; k++) {
			if (k == 0) {
				factor = sqrt(1.0 / block_size);
			}
			else {
				factor = sqrt(2.0 / block_size);
			}
			basis[n + k * block_size] = factor * cos((pi / block_size) * (n + 0.5) * k);
		}
	}
	return basis;
}

int Q[8 * 8] = { //given and optimized by hand
	   16,11,10,16,24,40,51,61,
	   12,12,14,19,26,58,60,55,
	   14,13,16,24,40,57,69,56,
	   14,17,22,29,51,87,80,62,
	   18,22,37,56,68,109,103,77,
	   24,35,55,64,81,104,113,92,
	   49,64,78,87,103,121,120,101,
	   72,92,95,98,112,100,103,99 };


float* dicing_window(float* image, int l, int k, int block_size) {
	float* image_as_block = new float[block_size * block_size];
	for (int i = 0; i < block_size; i++) {
		for (int j = 0; j < block_size; j++) {
			image_as_block[i + block_size * j] = image[l + i + ((k + j) * 256)];
		}
	}
	return image_as_block;
}
float* dicing_window1(float* image, int l, int k, int block_size) {
	float* image_as_block = new float[block_size * block_size];
	for (int i = 0; i < block_size; i++) {
		for (int j = 0; j < block_size; j++) {
			image_as_block[i + block_size * j] = image[l + i + ((k + j) * (256 - 8))];
		}
	}
	return image_as_block;
}
void insert_block(float* image, const float* block, int k, int l, int block_size) {

	for (int i = 0; i < block_size; i++) {
		for (int j = 0; j < block_size; j++) {
			//cout << block[i + block_size * j] << "   ";
			image[k * block_size + i + (l * block_size + j) * 256] = block[i + block_size * j];
		}
	}

}
void insert_block3(float* image, const float* block, int k, int l, int block_size) {

	for (int i = 0; i < block_size; i++) {
		for (int j = 0; j < block_size; j++) {
			//cout << block[i + block_size * j] << "   ";
			image[k * block_size + i + (l * block_size + j) * 256] = block[i *block_size+   j];
		}
	}

}
void insert_block1(float* image, const float* block, int k, int l, int block_size, int size) {

	for (int i = 0; i < block_size; i++) {
		for (int j = 0; j < block_size; j++) {
			image[k +(size/ block_size) * i + (l +( size/block_size) * j) * (256)] = block[i + block_size * j];
		}
	}

}
void insert_block2(float* image, const float* block, int k, int l, int block_size, int size) {

	for (int i = 0; i < block_size; i++) {
		for (int j = 0; j < block_size; j++) {
			image[k + block_size * i + (l + block_size * j) * (256)] = block[i + size/block_size * j];
		}
	}

}


void insert_dc_terms(float* image, const float* block, int k, int l, int block_size) {

	image[k + l * 32] = block[0];

}
float* delta_decoding(int size) {
	float* delta = new float[size * size];
	ifstream infile;
	infile.open("delta_encoding_dc.txt");
	float a;
	int l = 0;
	while (infile >> a)
	{
		delta[l] = a;
		l++;
	}
	float* dc_terms = new float[size * size];
	l = 0;
	for (int n = 0; n < size; ++n) {
		for (int k = 0; k < size; k++) {
			if (k == 0) {
				dc_terms[n + size * k] = delta[l];
				if (n == 0) {
					dc_terms[n + size * k] = delta[l];
				}
				else {
					dc_terms[n + size * k] = delta[l] + dc_terms[n - 1 + size * k];
				}

			}
			else {
				dc_terms[n + size * k] = delta[l] + dc_terms[n + size * (k - 1)];
			}
			l++;
		}

	}
	return dc_terms;

}
float* delta_decoding2(float*block,int size) {
	float* delta = new float[size * size];
	float* dc_terms = new float[size * size];
	for (int n = 0; n < size; ++n) {
		for (int k = 0; k < size; k++) {
			if (k == 0) {
				dc_terms[n + size * k] = block[n + size * k];
				if (n == 0) {
					dc_terms[n + size * k] = block[n + size * k];
				}
				else {
					dc_terms[n + size * k] = block[n + size * k] + dc_terms[n - 1 + size * k];
				}

			}
			else {
				dc_terms[n + size * k] = block[n + size * k] + dc_terms[n + size * (k - 1)];
			}
		}
	}
	return dc_terms;
}
float* delta_encoding(float* img_dc_terms, int size) {
	float* block = new float[size * size];
	ofstream myfile("delta_encoding_dc.txt");
	int l = 0;
	float* delta = new float[size * size];
	for (int n = 0; n < size; n++) {
		for (int k = 0; k < size; k++) {

			if (k == 0) {
				if (n == 0) {
					delta[l] = img_dc_terms[n + size * k];
					myfile << delta[l] << endl;
					block[n + k * size] = delta[l];
				}
				else {
					delta[l] = (img_dc_terms[n + size * k] - img_dc_terms[(n - 1) + size * k]);
					myfile << delta[l] << endl;
					block[n + k * size] = delta[l];
				}

			}
			else {
				delta[l] = (img_dc_terms[n + size * k] - img_dc_terms[n + size * (k - 1)]);
				myfile << delta[l] << endl;
				block[n + k * size] = delta[l];
			}
			l++;
		}
	}
	myfile.close();
	return block;
}
void run_length_encoding(float* img_ac_terms, int size, int block_size) {
	float* block = new float[block_size * block_size];
	int run_length = 0;
	int max_run_length = 0;
	ofstream myfile("run_length_encoding_ac_decompress.txt");
	long long int count_zeros = 0;
	for (int n = 0; n < size / block_size; n++) {
		for (int k = 0; k < size / block_size; k++) {
			block = dicing_window(img_ac_terms, block_size * n, block_size * k, block_size);
			count_zeros = 0;
			for (int x = 0; x < block_size; x++) {
				for (int y = 0; y < block_size; y++) {
					if ((block[x + block_size * y]) != 0) {
						if (count_zeros > 0) {
							myfile << 0 << endl;
							myfile << count_zeros << endl;
							myfile << block[x + block_size * y] << endl;
							run_length += 3;
							count_zeros = 0;
						}
						else {
							myfile << block[x + block_size * y] << endl;
							run_length++;
							count_zeros = 0;

						}

					}
					else {
						if (x == block_size - 1 && y == block_size - 1) {
							count_zeros += 1;
							myfile << 0 << endl;
							myfile << count_zeros << endl;
							run_length += 2;
							count_zeros = 0;


						}
						else {
							count_zeros += 1;

						}
					}
				}
			}
			count_zeros = 0;
			if (run_length > max_run_length && run_length < 1097) {
				max_run_length = run_length;
				cout << max_run_length << "    " << n<< k << endl;
			}
			run_length = 0;
			
		}
	}
	myfile.close();
	cout << "max run length =     " << max_run_length << endl;

}
float* run_length_decoding(const int size, const int block_size) {
	string* elemen = new string[size * size];
	ifstream infile;
	infile.open("run_length_encoding_ac_decompress.txt");
	string a;
	int l = 0;
	float* ac_term = new float[size * size];
	float* block = new float[block_size * block_size];
	while (getline(infile, a))
	{
		elemen[l] = a;
		l++;
	}
	infile.close();
	int i = 0;
	int index = 0;

	int idx = 0;
	int idy = 0;
	int count = 0;
	while (i < l) {
		if (count < block_size*block_size) {
			if (elemen[i] != "0") {
				block[index] = std::stoi(elemen[i]);
				count++;
				i++;
				index++;
			}
			else {
				for (int x = 0; x < std::stoi(elemen[i + 1]); x++) {
					
					block[index + x] = 0;
					count++;
					
				}
				index = index + std::stoi(elemen[i + 1]);
				i += 2;
			}
		}
		if (count == block_size * block_size) {
			insert_block3(ac_term, block, idx, idy, block_size);
			
			idy++;
			if (idy == (size / block_size)) {
				idy = 0; 
				idx++;
			}
			index = 0;
			count = 0;
		}
	}

	return ac_term;
}
void count_value_rle(int size) {
	string* elem = new string[size * size];
	ifstream ac;
	ac.open("run_length_encoding_ac_decompress.txt");
	string a;
	int l = 0;
	int total_symbol = 0;

	while (getline(ac, a))
	{
		elem[l] = a;
		l++;
	}
	total_symbol = l;
	map<string, int> dict;
	map<int, int>::iterator it;
	int i = 1097; 
	while (i < l) {
		++dict[elem[i]];
		i++;
	}
	int Histsize = dict.size();
	int totalelements = 0;
	for (auto iter = dict.begin(); iter != dict.end(); ++iter)
	{
		cout << "element: " << iter->first << "  frequency: " << iter->second << endl;
		totalelements = totalelements + (int)(iter->second);
	}
	string* elements = new string[Histsize];
	float* Probabilities = new float[Histsize];
	int n = 0;
	int intvalue;
	float Entropy = 0;
	for (auto iter = dict.begin(); iter != dict.end(); ++iter)
	{
		cout << "element: " << iter->first << "  frequency: " << iter->second;

		elements[n] = (iter->first);
		intvalue = (iter->second);
		Probabilities[n] = (float)((intvalue + 0.0) / (totalelements + 0.0));
		cout << " P(i) " << Probabilities[n] << endl;
		Entropy = Entropy - (Probabilities[n] * log2(Probabilities[n]));
		n++;
	}
	cout << "The total number of symbols is: " << totalelements << endl;
	cout << "The Number of symbols (Histogram) is: " << Histsize << endl;
	cout << "The Average bpsymbol is: " << Entropy << endl;
	cout << "The File size is: " << ((Entropy)*totalelements) / 8.0 << " bytes" << endl;
}
int inv_golomb(string binary_string) {
	int x = 0;
	int offset = 0;
	int N;
	int m = 0;
	int i = 0;
	bool stop = 1;
	int temp;
	while (binary_string[i] == '0') {
		m++;
		i++;
	}
	for (int k = 0; k < m; k++) {
		if (binary_string[2 * m - k] == '1') {
			temp = 1;

		}
		else {
			temp = 0;
		}

		offset += pow(2, k) * temp;
		//cout << offset << endl;
	}
	x = pow(2, m) - 1 + offset;


	if (x % 2 == 0) {
		x = -x / 2;
	}
	else {
		x = (x + 1) / 2;
	}
	return x;
}
string golomb(int x) {
	if (x > 0) {
		x = x * 2 - 1;
	}
	else {
		x = (-x) * 2;
	}
	string string_bit;
	int y = x + 1;
	int num = floor(log2(y)) + 1;
	int m = num - 1;
	int* array_of_bit = new int[num];
	while (m > 0) {
		string_bit += std::to_string(0);
		m--;
	}
	int i = 0;
	while (y > 0) {
		array_of_bit[i] = y % 2;
		y = y / 2;
		i++;
	}
	for (int j = i - 1; j >= 0; j--) {
		string_bit += std::to_string(array_of_bit[j]);
	}
	return string_bit;
}

void replace_dc_by_delta(float* img_ac ,int size, int block_size) {
	float* delta = new float[size * size];
	ifstream infile;
	infile.open("delta_encoding_dc.txt");
	float a;
	int l = 0;
	while (infile >> a)
	{
		delta[l] = a;
		l++;
	}
	float* dc_terms = new float[size/block_size * size / block_size];
	l = 0;
	for (int n = 0; n < size / block_size; ++n) {
		for (int k = 0; k < size / block_size; k++) {

			dc_terms[n + (size / block_size) * k] = delta[l];
		}
	}
	for (int n = 0; n < size / block_size; n ++) {
		for (int k = 0; k < size / block_size; k ++) {
			img_ac[(n * size / block_size) + size * (k * size / block_size)] = dc_terms[n + (size / block_size )* k];
		}
	}
}
void compress(float* img, int img_size) {
	const int block_size = 8;
	const float* basis = get_basis(block_size);
	float* img_dc_term = new float[(img_size / block_size) * (img_size / block_size)];
	float* img_ac_term = new float[(img_size) * (img_size)];
	for (int x = 0; x < (img_size / block_size); x++) {
		for (int y = 0; y < (img_size / block_size); y++) {
			float* block = dicing_window(img, x * block_size, y * block_size, block_size);
			float* block_approx = encode(block, basis, Q, block_size);
			insert_block1(img_ac_term, block_approx, x, y, block_size, img_size);
			
		}
	}
	Store(img_ac_term, "img_ac_term.raw", img_size);
	img_dc_term = dicing_window(img_ac_term, 0, 0, img_size / block_size);
	float*img_dc_term_delta = delta_encoding(img_dc_term, (img_size / block_size));
	insert_block(img_ac_term, img_dc_term_delta, 0, 0, img_size / block_size);
	Store(img_ac_term, "img_ac_term_dc_delta.raw", img_size);
	run_length_encoding(img_ac_term, img_size,(img_size/ block_size));
	string* elemen = new string[img_size * img_size];
	ifstream infile;
	infile.open("run_length_encoding_ac_decompress.txt");
	string a;
	int i = 0;
	int l = 0;
	int ac_term = 0;
	ofstream myfile("compressed.txt");
	while (getline(infile, a))
	{
		elemen[l] = a;
		l++;
	}
	infile.close();
	while (i < l) {
		ac_term = std::stoi(elemen[i]);
		i++;
		myfile << golomb(ac_term) << endl;
	}

	myfile.close();
	

}
float* decompress(int img_size,int  block_size) {
	ifstream infile;
	string a;
	string* elemen = new string[img_size * img_size];
	
	int i = 0;
	int l = 0;
	infile.open("compressed.txt");
	while (getline(infile, a))
	{
		elemen[l] = a;
		l++;
	}
	infile.close();
	ofstream myfile("run_length_encoding_ac_decompress.txt");
	while (i < l) {
		myfile << inv_golomb(elemen[i]) << endl;
		i++;
	}
	myfile.close();
	
	float* img_ac_term_reconstruct = run_length_decoding(img_size, img_size/block_size);
	Store(img_ac_term_reconstruct, "img_ac_term_reconstruct_without_delta_decoding.raw", img_size);
	float* dc_block = dicing_window(img_ac_term_reconstruct, 0, 0, img_size / block_size);
	Store(dc_block, "dc_block.raw", img_size / block_size);
	float* dc_block_decoded = delta_decoding2(dc_block, img_size / block_size);
	insert_block(img_ac_term_reconstruct, dc_block_decoded, 0, 0, img_size / block_size);
	Store(img_ac_term_reconstruct, "img_ac_term_reconstruct_with_delta_decoding.raw", img_size);

	float* layout_for_idct = new float[img_size * img_size];

	for (int x = 0; x < block_size; x++) {
		for (int y = 0; y <  block_size; y++) {
			float* block = dicing_window(img_ac_term_reconstruct, x * img_size / block_size, y * img_size / block_size, img_size/block_size);
			insert_block1(layout_for_idct, block, x, y, img_size/block_size, img_size);
		}
	}
	Store(layout_for_idct, "layout_for_idct.raw", img_size);
	float* img_retrivial = new float[img_size * img_size];
	for (int x = 0; x < (img_size / block_size); x++) {
		for (int y = 0; y < (img_size / block_size); y++) {
			float* block = dicing_window(layout_for_idct, x * block_size, y * block_size, block_size);
			float* block_approx = decode(block,get_basis(block_size), Q, block_size);
			insert_block(img_retrivial, block_approx, x, y, block_size);

		}
	}
	return img_retrivial;
	

}


int main(){
	float* lena = load("lena_256x256.raw");
	int img_size = 256;
	int block_size = 8;
	compress(lena, img_size);
	float* img_decompressed = decompress(img_size, block_size);
	Store(img_decompressed, "img_decompressed.raw", img_size);

	float mse; 
	float psnr; 

	mse = MSE(lena, img_decompressed, 256);
	psnr = PSNR(255, lena, img_decompressed, 256);
	
	std::cout << "mse = " << mse << std::endl;
	std::cout << "psnr = " << psnr << std::endl;
	count_value_rle(256);
	int i = -20;
	for (i; i < 20; i++) {
		cout << "number:  " << i << "     codeword golomb:   " << (golomb(i)) << "    result -->  " << inv_golomb(golomb(i)) << endl;
	}
}
