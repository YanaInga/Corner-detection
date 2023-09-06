#include <iostream>
#include <omp.h>
#include <string>
#include <Windows.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <tbb/tbb.h>
using namespace std;


void BMPWrite(RGBQUAD**&, BITMAPFILEHEADER&, BITMAPINFOHEADER&, const char*);
void BMPRead(RGBQUAD**&, BITMAPFILEHEADER&, BITMAPINFOHEADER&, const char*);


void BMPRead(RGBQUAD**& rgb, BITMAPFILEHEADER& header, \
    BITMAPINFOHEADER& bmiHeader, const char* fin)
{
    // Открываем файл
    std::ifstream InFile(fin, std::ios::binary);
    // Считываем заголовок файла
    InFile.read((char*)(&header), sizeof(BITMAPFILEHEADER));
    // Считываем заголовочную часть изображения
    InFile.read((char*)(&bmiHeader), sizeof(BITMAPINFOHEADER));
    // Выделяем память под массив RGB хранящий структуры RGBQUAD
    rgb = new RGBQUAD * [bmiHeader.biHeight];
    for (int i = 0; i < bmiHeader.biHeight; i++)
    {
        rgb[i] = new RGBQUAD[bmiHeader.biWidth];
    }
    // Считываем данные изображения в массив структур RGB 
    for (int i = 0; i < bmiHeader.biHeight; i++)
    {
        for (int j = 0; j < bmiHeader.biWidth; j++)
        {
            InFile.read((char*)(&rgb[i][j]), 3); // .rgbBlue .rgbGreen .rgbRed;
        }
    }
    // Закрываем файл
    InFile.close();

}

void BMPWrite(RGBQUAD**& rgb, BITMAPFILEHEADER& header, \
    BITMAPINFOHEADER& bmiHeader, const char* fout)
{
    // Открываем файл для записи изображения в формат BMP
    std::ofstream OutFile(fout, std::ios::binary);
    //// Записываем заголовок файла
    OutFile.write((char*)(&header), sizeof(BITMAPFILEHEADER));
    //// Записываем заголовочную часть изображения
    OutFile.write((char*)(&bmiHeader), sizeof(BITMAPINFOHEADER));
    // Записываем данные изображения из массив структур RGB в файл 
    for (int i = 0; i < bmiHeader.biHeight; i++)
        for (int j = 0; j < bmiHeader.biWidth; j++)
        {
            OutFile.write((char*)&(rgb[i][j]), 3);// .rgbBlue .rgbGreen .rgbRed;
        }
    // закрываем файл
    OutFile.close();
}
class HarrisPoints
{
public:
    double R;
    int x;
    int y;
    HarrisPoints(double _R, int _x, int _y) : R(_R), x(_x), y(_y) {}
};

void findcoeff(double**& MATR_coef, int RH, int RW, int ksize)
{
	double sum = 0;
	for (int y = -RH; y <= RH; y++)
		for (int x = -RW; x <= RW; x++) {
			MATR_coef[y + RH][x + RW] = (1 / (2 * 3.14 * ksize * ksize)) * exp(-1 * (x * x + y * y) / (ksize * ksize * 2));
			sum += MATR_coef[y + RH][x + RW];
		}
	for (int y = -RH; y <= RH; y++)
		for (int x = -RW; x <= RW; x++) {
			MATR_coef[y + RH][x + RW] /= sum;
		}
}
void gaussfilter(int ksize, int height, int width, double** rgb1, double**& rgb2)
{
	int size = ksize * ksize;
	int rh = (ksize - 1) / 2, rw = (ksize - 1) / 2;
	double** Matr_coef = new double* [ksize];
	for (int i = 0; i < ksize; i++)
		Matr_coef[i] = new double[ksize];
	findcoeff(Matr_coef, rh, rw, ksize);
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			double LinF_Value = 0.0;
			for (int dy = -rh; dy <= rh; dy++)
			{
				int ky = y + dy;
				if (ky < 0) ky = 0;
				if (ky > height - 1) ky = height - 1;
				for (int dx = -rw; dx <= rw; dx++)
				{
					int kx = x + dx;
					if (kx < 0) kx = 0;
					if (kx > width - 1) kx = width - 1;
					LinF_Value += rgb1[ky][kx] * Matr_coef[dy + rh][dx + rw];
				}
			}
			if (LinF_Value < 0) LinF_Value = 0;
			if (LinF_Value > 255) LinF_Value = 255;
			rgb2[y][x] = (int)LinF_Value;
		}
	}
}
void gaussfilterp(int ksize, int height, int width, double** rgb1, double**& rgb2)
{
	int size = ksize * ksize;
	int rh = (ksize - 1) / 2, rw = (ksize - 1) / 2;
	double** Matr_coef = new double* [ksize];
	for (int i = 0; i < ksize; i++)
		Matr_coef[i] = new double[ksize];
	findcoeff(Matr_coef, rh, rw, ksize);
#pragma omp parallel for
	for (int i = 0; i < ksize; i++)
		Matr_coef[i] = new double[ksize];
	findcoeff(Matr_coef, rh, rw, ksize);
#pragma omp parallel for schedule(dynamic, height/12)
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			double LinF_Valuer = 0.0;
			for (int dy = -rh; dy <= rh; dy++)
			{
				int ky = y + dy;
				if (ky < 0) ky = 0;
				if (ky > height - 1) ky = height - 1;
				for (int dx = -rw; dx <= rw; dx++)
				{
					int kx = x + dx;
					if (kx < 0) kx = 0;
					if (kx > width - 1) kx = width - 1;
					LinF_Valuer += rgb1[ky][kx] * Matr_coef[dy + rh][dx + rw];

				}
			}
			if (LinF_Valuer < 0) LinF_Valuer = 0;
			if (LinF_Valuer > 255) LinF_Valuer = 255;
			rgb2[y][x] = (int)LinF_Valuer;
		}
	}

}
void gaussfilteri(int ksize, int height, int width, double** rgb1, double**& rgb2)
{
	int size = ksize * ksize;
	int rh = (ksize - 1) / 2, rw = (ksize - 1) / 2;
	double** Matr_coef = new double* [ksize];
	for (int i = 0; i < ksize; i++)
		Matr_coef[i] = new double[ksize];
	findcoeff(Matr_coef, rh, rw, ksize);
	tbb::parallel_for(tbb::blocked_range2d<int>(0, height, 0, width), [&]
	(tbb::blocked_range2d<int> r)
		{
			for (int y = r.rows().begin(); y < r.rows().end(); y++)
				for (int x = r.cols().begin(); x < r.cols().end(); x++)
				{
					double LinF_Valuer = 0.0;
					double LinF_Valueg = 0.0;
					double LinF_Valueb = 0.0;
					for (int dy = -rh; dy <= rh; dy++)
					{
						int ky = y + dy;
						if (ky < 0) ky = 0;
						if (ky > height - 1) ky = height - 1;
						for (int dx = -rw; dx <= rw; dx++)
						{
							int kx = x + dx;
							if (kx < 0) kx = 0;
							if (kx > width - 1) kx = width - 1;
							LinF_Valuer += rgb1[ky][kx] * Matr_coef[dy + rh][dx + rw];
						}
					}
					if (LinF_Valuer < 0) LinF_Valuer = 0;
					if (LinF_Valuer > 255) LinF_Valuer = 255;
					rgb2[y][x] = (int)LinF_Valuer;
				}
		});
}
template <class T>
class MaxCalc
{
private:
	T** ArrayA;
public:
	T MaxValue;
	T MinValue;
	void operator()(const tbb::blocked_range2d<size_t>& r)
	{
		T** a = ArrayA;
		T val;
		for (int i = r.rows().begin(); i != r.rows().end(); ++i)
			for (int j = r.cols().begin(); j != r.cols().end(); ++j)
			{
				val = a[i][j];
				if (val < MinValue)MinValue = val;
				if (val > MaxValue)MaxValue = val;
			}
	}
	MaxCalc<T>(MaxCalc& x, tbb::split) :
		ArrayA(x.ArrayA), MinValue(DBL_MAX), MaxValue(-DBL_MAX) {}
	void join(const MaxCalc& y)
	{
		if (y.MinValue < MinValue) MinValue = y.MinValue;
		if (y.MaxValue > MaxValue) MaxValue = y.MaxValue;
	}
	MaxCalc<T>(T** A) :
		ArrayA(A),
		MinValue(DBL_MAX),
		MaxValue(-DBL_MAX)// Инициализация MaxValue
	{}
};
