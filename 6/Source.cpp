#include <iostream>
#include <omp.h>
#include <string>
#include <Windows.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include "BMPFileRW.h"
#include <fstream>
#include <tbb/tbb.h>
using namespace std;

double ArmV(double* time) //Среднее арифметическое значение
{
	double amvTime = 0;
	for (int i = 0; i < 10; i++)
	{
		amvTime += time[i];
	}
	return (double)amvTime / 10;
}
// Функция рассчета среднеарифметического значения в доверительном интервале
double AvgTrustedInterval(double& amvTime, double* time)
{
	double sd = 0, newAVg = 0;
	int newiter = 0;
	for (int i = 0; i < 10; i++)
	{
		sd += (time[i] - amvTime) * (time[i] - amvTime);
	}
	sd /= 9.0;
	sd = sqrt(sd);
	for (int i = 0; i < 10; i++)
	{
		if (amvTime - sd <= time[i] && time[i] <= amvTime + sd)
		{
			newAVg += time[i];
			newiter++;
		}
	}
	if (newiter == 0) newiter = 1;
	return newAVg / newiter;
}
void bright(int height, int width, RGBQUAD** rgb, double** I, double** Inew)
{
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
		{
			I[y][x] = rgb[y][x].rgbRed * 0.299 + rgb[y][x].rgbGreen * 0.587 + rgb[y][x].rgbBlue * 0.114;
		}
	gaussfilter(5, height, width, I, Inew);
}
void gradient(int height, int width, double** I, double** DiffX, double** DiffY, double** DiffXY)
{
	int A[3][3] = { { 1, 0, -1 },
					{ 1, 0, -1 },
					{ 1, 0, -1 } };
	double** hor = new double* [height];
	double** ver = new double* [height];
	for (int i = 0; i < height; i++)
	{
		hor[i] = new double[width];
		ver[i] = new double[width];
	}
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
		{
			hor[y][x] = 0;
			ver[y][x] = 0;
			for (int k = 0; k <= 2; k++)
			{
				int dk = y + k - 1;
				if (dk < 0) dk = 0;
				if (dk > height - 1) dk = height - 1;
				for (int l = 0; l <= 2; l++)
				{
					int dl = x + l - 1;
					if (dl < 0) dl = 0;
					if (dl > width - 1) dl = width - 1;
					hor[y][x] += A[k][l] * I[dk][dl];
					ver[y][x] += A[l][k] * I[dk][dl];
				}
			}
			if (hor[y][x] < 0)
				hor[y][x] = 0;
			if (hor[y][x] > 255)
				hor[y][x] = 255;
			if (ver[y][x] < 0)
				ver[y][x] = 0;
			if (ver[y][x] > 255)
				ver[y][x] = 255;
		}
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
		{
			DiffX[y][x] = 0;
			DiffY[y][x] = 0;
			DiffXY[y][x] = 0;
			for (int k = 0; k <= 2; k++)
			{
				int dk = y + k - 1;
				if (dk < 0) dk = 0;
				if (dk > height - 1) dk = height - 1;
				for (int l = 0; l <= 2; l++)
				{
					int dl = x + l - 1;
					if (dl < 0) dl = 0;
					if (dl > width - 1) dl = width - 1;
					DiffX[y][x] += A[k][l] * hor[dk][dl];
					DiffY[y][x] += A[l][k] * ver[dk][dl];
					DiffXY[y][x] += A[l][k] * hor[dk][dl];
				}
			}
			if (DiffX[y][x] < 0)
				DiffX[y][x] = 0;
			if (DiffX[y][x] > 255)
				DiffX[y][x] = 255;
			if (DiffY[y][x] < 0)
				DiffY[y][x] = 0;
			if (DiffY[y][x] > 255)
				DiffY[y][x] = 255;
			if (DiffXY[y][x] < 0)
				DiffXY[y][x] = 0;
			if (DiffXY[y][x] > 255)
				DiffXY[y][x] = 255;
		}
}
void funcR(int height, int width, double** I, double** H, double** R, double** DiffX, double** DiffY, double** DiffXY, double threshold)
{
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
		{
			R[y][x] = (DiffX[y][x] * DiffY[y][x] - DiffXY[y][x] * DiffXY[y][x]) - (0.04 * (DiffX[y][x] + DiffY[y][x]) * (DiffX[y][x] + DiffY[y][x]));
		}
	double min = R[0][0], max = R[0][0];
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
		{
			if (R[y][x] < min) min = R[y][x];
			if (R[y][x] > max) max = R[y][x];
		}
	threshold = ((max - min) * threshold) / 100;
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
		{
			if (R[y][x] > threshold)
				H[y][x] = R[y][x];
			else H[y][x] = 0;
		}
}
void findpoints(int height, int width, double** H, vector<HarrisPoints>& HP, int rang)
{
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
		{
			if (H[y][x] > 0)
			{
				HarrisPoints i (H[y][x], x, y);
				HP.push_back(i);
			}
		}
	for (int i = 0; i < HP.size(); i++)
		for (int j = HP.size()-1; j > i; j--)
		{
			if (HP[j-1].R < HP[j].R)
			{
				swap(HP[j-1], HP[j]);
			}
		}
	for (int i = 0; i < HP.size(); i++)
		for (int j = 0; j < HP.size(); j++)
		{
			if((HP[i].x - rang <= HP[j].x) && (HP[i].x + rang >= HP[j].x)
				&& (HP[i].y - rang <= HP[j].y) && (HP[i].y + rang >= HP[j].y))
			{
				HP.erase(HP.begin() + j);
				H[HP[j].y][HP[j].x] = 0; 
			}
		}
}
void drawpoints(int height, int width,vector<HarrisPoints>& HP, RGBQUAD**& rgbres)
{
	for (int i = 0; i < HP.size(); i++)
	{
		for (int y = -1; y < 2; y++)
		{
			for (int x = -1; x < 2; x++)
			{
				int dy = HP[i].y + y;
				int dx = HP[i].x + x;

				if (dy >= 0 && dy <= height - 1 && dx >= 0 && dx <= width - 1)
				{
					rgbres[dy][dx].rgbRed = 255;
					rgbres[dy][dx].rgbGreen = 0;
					rgbres[dy][dx].rgbBlue = 0;
				}
			}
		}
	}
}
void brightp(int height, int width, RGBQUAD** rgb, double** I, double** Inew)
{
#pragma omp parallel for
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
		{
			I[y][x] = rgb[y][x].rgbRed * 0.299 + rgb[y][x].rgbGreen * 0.587 + rgb[y][x].rgbBlue * 0.114;
		}
	gaussfilterp(5, height, width, I, Inew);
}
void gradientp(int height, int width, double** I, double** DiffX, double** DiffY, double** DiffXY)
{
		int A[3][3] = { { 1, 0, -1 },
						{ 1, 0, -1 },
						{ 1, 0, -1 } };
		double** hor = new double* [height];
		double** ver = new double* [height];
#pragma omp parallel for schedule(dynamic, height/12)
		for (int i = 0; i < height; i++)
		{
			hor[i] = new double[width];
			ver[i] = new double[width];
		}
#pragma omp parallel for schedule(dynamic, height/12)
		for (int y = 0; y < height; y++)
			for (int x = 0; x < width; x++)
			{
				hor[y][x] = 0;
				ver[y][x] = 0;
				for (int k = 0; k <= 2; k++)
				{
					int dk = y + k - 1;
					if (dk < 0) dk = 0;
					if (dk > height - 1) dk = height - 1;
					for (int l = 0; l <= 2; l++)
					{
						int dl = x + l - 1;
						if (dl < 0) dl = 0;
						if (dl > width - 1) dl = width - 1;
						hor[y][x] += A[k][l] * I[dk][dl];
						ver[y][x] += A[l][k] * I[dk][dl];
					}
				}
				if (hor[y][x] < 0)
					hor[y][x] = 0;
				if (hor[y][x] > 255)
					hor[y][x] = 255;
				if (ver[y][x] < 0)
					ver[y][x] = 0;
				if (ver[y][x] > 255)
					ver[y][x] = 255;
			}
#pragma omp parallel for schedule(dynamic, height/12)
		for (int y = 0; y < height; y++)
			for (int x = 0; x < width; x++)
			{
				DiffX[y][x] = 0;
				DiffY[y][x] = 0;
				DiffXY[y][x] = 0;
				for (int k = 0; k <= 2; k++)
				{
					int dk = y + k - 1;
					if (dk < 0) dk = 0;
					if (dk > height - 1) dk = height - 1;
					for (int l = 0; l <= 2; l++)
					{
						int dl = x + l - 1;
						if (dl < 0) dl = 0;
						if (dl > width - 1) dl = width - 1;
						DiffX[y][x] += A[k][l] * hor[dk][dl];
						DiffY[y][x] += A[l][k] * ver[dk][dl];
						DiffXY[y][x] += A[l][k] * hor[dk][dl];
					}
				}
				if (DiffX[y][x] < 0)
					DiffX[y][x] = 0;
				if (DiffX[y][x] > 255)
					DiffX[y][x] = 255;
				if (DiffY[y][x] < 0)
					DiffY[y][x] = 0;
				if (DiffY[y][x] > 255)
					DiffY[y][x] = 255;
				if (DiffXY[y][x] < 0)
					DiffXY[y][x] = 0;
				if (DiffXY[y][x] > 255)
					DiffXY[y][x] = 255;
			}
}
void funcRp(int height, int width, double** I, double** H, double** R, double** DiffX, double** DiffY, double** DiffXY, double threshold)
{
#pragma omp parallel for schedule(dynamic, height/12)
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
		{
			R[y][x] = (DiffX[y][x] * DiffY[y][x] - DiffXY[y][x] * DiffXY[y][x]) - (0.04 * (DiffX[y][x] + DiffY[y][x]) * (DiffX[y][x] + DiffY[y][x]));
		}
	double min = R[0][0], max = R[0][0];
#pragma omp parallel for schedule(dynamic, height/12) reduction (min:min) reduction(max:max)
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
		{
			if (R[y][x] < min) min = R[y][x];
			if (R[y][x] > max) max = R[y][x];
		}
	threshold = ((max - min) * threshold) / 100;
#pragma omp parallel for shared(threshold)
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
		{
			if (R[y][x] > threshold)
				H[y][x] = R[y][x];
			else H[y][x] = 0;
		}
}
void findpointsp(int height, int width, double** H, vector<HarrisPoints>& HP, int rang)
{
#pragma omp parallel for schedule(dynamic, HP.size()/12)
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
		{
			if (H[y][x] > 0)
			{
				HarrisPoints i(H[y][x], x, y);
#pragma omp critical
				HP.push_back(i);
			}
		}
	for (int i = 0; i < HP.size(); i++)
		for (int j = HP.size() - 1; j > i; j--)
		{
			if (HP[j - 1].R < HP[j].R)
			{
				swap(HP[j - 1], HP[j]);
			}
		}
	for (int i = 0; i < HP.size(); i++)
		for (int j = 0; j < HP.size(); j++)
		{
			if ((HP[i].x - rang <= HP[j].x) && (HP[i].x + rang >= HP[j].x)
				&& (HP[i].y - rang <= HP[j].y) && (HP[i].y + rang >= HP[j].y))
			{
				HP.erase(HP.begin() + j);
				H[HP[j].y][HP[j].x] = 0;
			}
		}
}
void drawpointsp(int height, int width, vector<HarrisPoints>& HP, RGBQUAD**& rgbres)
{
#pragma omp parallel for schedule(dynamic, HP.size()/12)
	for (int i = 0; i < HP.size(); i++)
	{
		for (int y = -1; y < 2; y++)
		{
			for (int x = -1; x < 2; x++)
			{
				int dy = HP[i].y + y;
				int dx = HP[i].x + x;

				if (dy >= 0 && dy <= height - 1 && dx >= 0 && dx <= width - 1)
				{
					rgbres[dy][dx].rgbRed = 255;
					rgbres[dy][dx].rgbGreen = 0;
					rgbres[dy][dx].rgbBlue = 0;
				}
			}
		}
	}
}
void brighti(int height, int width, RGBQUAD** rgb, double** I, double** Inew)
{
	tbb::parallel_for(tbb::blocked_range2d<int>(0, height, 0, width), [&]
	(tbb::blocked_range2d<int> r)
		{
			for (int y = r.rows().begin(); y < r.rows().end(); y++)
				for (int x = r.cols().begin(); x < r.cols().end(); x++)
				{
					I[y][x] = rgb[y][x].rgbRed * 0.299 + rgb[y][x].rgbGreen * 0.587 + rgb[y][x].rgbBlue * 0.114;
				}
		});
	gaussfilteri(5, height, width, I, Inew);
}
void gradienti(int height, int width, double** I, double** DiffX, double** DiffY, double** DiffXY)
{
	int A[3][3] = { { 1, 0, -1 },
					{ 1, 0, -1 },
					{ 1, 0, -1 } };
	double** hor = new double* [height];
	double** ver = new double* [height];
	tbb::parallel_for(tbb::blocked_range<int>(0, height), [&]
	(tbb::blocked_range<int> r)
		{for (int i = r.begin(); i < r.end(); i++)
	{
		hor[i] = new double[width];
		ver[i] = new double[width];
	}});
	tbb::parallel_for(tbb::blocked_range2d<int>(0, height, 0, width), [&]
	(tbb::blocked_range2d<int> r) {
			for (int y = r.rows().begin(); y < r.rows().end(); y++)
				for (int x = r.cols().begin(); x < r.cols().end(); x++)
				{
					hor[y][x] = 0;
					ver[y][x] = 0;
					for (int k = 0; k <= 2; k++)
					{
						int dk = y + k - 1;
						if (dk < 0) dk = 0;
						if (dk > height - 1) dk = height - 1;
						for (int l = 0; l <= 2; l++)
						{
							int dl = x + l - 1;
							if (dl < 0) dl = 0;
							if (dl > width - 1) dl = width - 1;
							hor[y][x] += A[k][l] * I[dk][dl];
							ver[y][x] += A[l][k] * I[dk][dl];
						}
					}
					if (hor[y][x] < 0)
						hor[y][x] = 0;
					if (hor[y][x] > 255)
						hor[y][x] = 255;
					if (ver[y][x] < 0)
						ver[y][x] = 0;
					if (ver[y][x] > 255)
						ver[y][x] = 255;
				}
		});
	tbb::parallel_for(tbb::blocked_range2d<int>(0, height, 0, width), [&]
	(tbb::blocked_range2d<int> r) {
			for (int y = r.rows().begin(); y < r.rows().end(); y++)
				for (int x = r.cols().begin(); x < r.cols().end(); x++)
				{
					DiffX[y][x] = 0;
					DiffY[y][x] = 0;
					DiffXY[y][x] = 0;
					for (int k = 0; k <= 2; k++)
					{
						int dk = y + k - 1;
						if (dk < 0) dk = 0;
						if (dk > height - 1) dk = height - 1;
						for (int l = 0; l <= 2; l++)
						{
							int dl = x + l - 1;
							if (dl < 0) dl = 0;
							if (dl > width - 1) dl = width - 1;
							DiffX[y][x] += A[k][l] * hor[dk][dl];
							DiffY[y][x] += A[l][k] * ver[dk][dl];
							DiffXY[y][x] += A[l][k] * hor[dk][dl];
						}
					}
					if (DiffX[y][x] < 0)
						DiffX[y][x] = 0;
					if (DiffX[y][x] > 255)
						DiffX[y][x] = 255;
					if (DiffY[y][x] < 0)
						DiffY[y][x] = 0;
					if (DiffY[y][x] > 255)
						DiffY[y][x] = 255;
					if (DiffXY[y][x] < 0)
						DiffXY[y][x] = 0;
					if (DiffXY[y][x] > 255)
						DiffXY[y][x] = 255;
				}});
}
void funcRi(int height, int width, double** I, double** H, double** R, double** DiffX, double** DiffY, double** DiffXY, double threshold)
{
	tbb::parallel_for(tbb::blocked_range2d<int>(0, height, 0, width), [&]
	(tbb::blocked_range2d<int> r) {
			for (int y = r.rows().begin(); y < r.rows().end(); y++)
				for (int x = r.cols().begin(); x < r.cols().end(); x++)
				{
					R[y][x] = (DiffX[y][x] * DiffY[y][x] - DiffXY[y][x] * DiffXY[y][x]) - (0.04 * (DiffX[y][x] + DiffY[y][x]) * (DiffX[y][x] + DiffY[y][x]));
				}});
	MaxCalc<double> MMC(R);
	tbb::parallel_reduce(tbb::blocked_range2d<size_t>(0, height, 0, width), MMC);
	double max = MMC.MaxValue;
	double min = MMC.MinValue;
	threshold = ((max - min) * threshold) / 100;
	tbb::parallel_for(tbb::blocked_range2d<int>(0, height, 0, width), [&]
	(tbb::blocked_range2d<int> r) {
			for (int y = r.rows().begin(); y < r.rows().end(); y++)
				for (int x = r.cols().begin(); x < r.cols().end(); x++)
				{
					if (R[y][x] > threshold)
						H[y][x] = R[y][x];
					else H[y][x] = 0;
				}});
}
void findpointsi(int height, int width, double** H, vector<HarrisPoints>& HP, int rang)
{
	mutex m;
	tbb::parallel_for(tbb::blocked_range2d<int>(0, height, 0, width), [&]
	(tbb::blocked_range2d<int> r) {
			for (int y = r.rows().begin(); y < r.rows().end(); y++)
				for (int x = r.cols().begin(); x < r.cols().end(); x++)
				{
					if (H[y][x] > 0)
					{
						HarrisPoints i(H[y][x], x, y);
						m.lock();
						HP.push_back(i);
						m.unlock();
					}
				}});
	for (int i = 0; i < HP.size(); i++)
		for (int j = HP.size() - 1; j > i; j--)
		{
			if (HP[j - 1].R < HP[j].R)
			{
				swap(HP[j - 1], HP[j]);
			}
		}
	for (int i = 0; i < HP.size(); i++)
		for (int j = 0; j < HP.size(); j++)
		{
			if ((HP[i].x - rang <= HP[j].x) && (HP[i].x + rang >= HP[j].x)
				&& (HP[i].y - rang <= HP[j].y) && (HP[i].y + rang >= HP[j].y))
			{
				HP.erase(HP.begin() + j);
				H[HP[j].y][HP[j].x] = 0;
			}
		}
}
void drawpointsi(int height, int width, vector<HarrisPoints>& HP, RGBQUAD**& rgbres)
{
	tbb::parallel_for(tbb::blocked_range<int>(0, HP.size()), [&]
	(tbb::blocked_range<int> r)
		{
			for (int i = r.begin(); i < r.end(); i++)
			{
				for (int y = -1; y < 2; y++)
				{
					for (int x = -1; x < 2; x++)
					{
						int dy = HP[i].y + y;
						int dx = HP[i].x + x;

						if (dy >= 0 && dy <= height - 1 && dx >= 0 && dx <= width - 1)
						{
							rgbres[dy][dx].rgbRed = 255;
							rgbres[dy][dx].rgbGreen = 0;
							rgbres[dy][dx].rgbBlue = 0;
						}
					}
				}

			}
		});
}
void formpictures(RGBQUAD**& rgb, int width, int height, BITMAPFILEHEADER header, BITMAPINFOHEADER bmiHeader)
{
	double* time = new double[10];
	double* time1 = new double[10];
	double* time2 = new double[10];
	double* time3 = new double[10];
	double* time4 = new double[10];
	double* time5 = new double[10];
	double** I = new double* [height];
	double** Inew = new double* [height];
	double** DiffX = new double* [height];
	double** DiffY = new double* [height];
	double** DiffXY = new double* [height];
	double** H = new double* [height];
	double** R = new double* [height];
	vector<HarrisPoints> HP;
	for (int i = 0; i < height; i++)
	{
		I[i] = new double[width];
		Inew[i] = new double[width];
		DiffX[i] = new double[width];
		DiffY[i] = new double[width];
		DiffXY[i] = new double[width];
		H[i] = new double[width];
		R[i] = new double[width];
	}
	double amvTime;
	//последовательно
	ofstream b, g, fR, fp, dp;
	b.open("c:\\test2\\bright.txt", ios::app);
	g.open("c:\\test2\\g.txt", ios::app);
	fR.open("c:\\test2\\fr.txt", ios::app);
	fp.open("c:\\test2\\fp.txt", ios::app);
	dp.open("c:\\test2\\dp.txt", ios::app);
	for (int i = 0; i < 10; i++)
	{
		double t_start = omp_get_wtime();
		//Формирование массива яркости изображения
		double t = omp_get_wtime();
		bright(height, width, rgb, I, Inew);
		double te = omp_get_wtime();
		time1[i] = (te - t) * 1000;
		//Вычисление первого и второго градиента яркости по горизонтали и вертикали
		t = omp_get_wtime();
		gradient(height, width, Inew, DiffX, DiffY, DiffXY);
		te = omp_get_wtime();
		time2[i] = (te - t) * 1000;
		//Вычисление значения функции отклика угла R
		t = omp_get_wtime();
		funcR(height, width, Inew, H, R, DiffX, DiffY, DiffXY, 20.0);
		te = omp_get_wtime();
		time3[i] = (te - t) * 1000;
		//Поиск особых точек для визуализации
		t = omp_get_wtime();
		findpoints(height, width, H, HP, 20);
		te = omp_get_wtime();
		time4[i] = (te - t) * 1000;
		//Визуализация особых точек
		t = omp_get_wtime();
		drawpoints(height, width, HP, rgb);
		te = omp_get_wtime();
		time5[i] = (te - t) * 1000;
		double t_end = omp_get_wtime();
		time[i] = (t_end - t_start) * 1000;		
	}
	amvTime = ArmV(time1);
	amvTime = AvgTrustedInterval(amvTime, time1);
	b << amvTime << endl;
	amvTime = ArmV(time2);
	amvTime = AvgTrustedInterval(amvTime, time2);
	g << amvTime << endl;
	amvTime = ArmV(time3);
	amvTime = AvgTrustedInterval(amvTime, time3);
	fR << amvTime << endl;
	amvTime = ArmV(time4);
	amvTime = AvgTrustedInterval(amvTime, time4);
	fp << amvTime << endl;
	amvTime = ArmV(time5);
	amvTime = AvgTrustedInterval(amvTime, time5);
	dp << amvTime << endl;
	amvTime = ArmV(time);
	amvTime = AvgTrustedInterval(amvTime, time);
	//cout << "Последовательно ";
	cout << amvTime << endl;
	string name = "c:\\test2\\posl" + to_string(width) + "x" + to_string(height) + ".bmp";
	BMPWrite(rgb, header, bmiHeader, name.c_str());
	//параллельно
	int count = 2;
	while (count < 5)
	{
		omp_set_num_threads(count);
		for (int i = 0; i < 10; i++)
		{
			double t_start = omp_get_wtime();
			//Формирование массива яркости изображения
			double t = omp_get_wtime();
			brightp(height, width, rgb, I, Inew);
			double te = omp_get_wtime();
			time1[i] = (te - t) * 1000;
			//Вычисление первого и второго градиента яркости по горизонтали и вертикали
			t = omp_get_wtime();
			gradientp(height, width, Inew, DiffX, DiffY, DiffXY);
			te = omp_get_wtime();
			time2[i] = (te - t) * 1000;
			//Вычисление значения функции отклика угла R
			t = omp_get_wtime();
			funcRp(height, width, Inew, H, R, DiffX, DiffY, DiffXY, 20.0);
			te = omp_get_wtime();
			time3[i] = (te - t) * 1000;
			//Поиск особых точек для визуализации
			t = omp_get_wtime();
			findpointsp(height, width, H, HP, 20);
			te = omp_get_wtime();
			time4[i] = (te - t) * 1000;
			//Визуализация особых точек
			t = omp_get_wtime();
			drawpointsp(height, width, HP, rgb);
			te = omp_get_wtime();
			time5[i] = (te - t) * 1000;
			double t_end = omp_get_wtime();
			time[i] = (t_end - t_start) * 1000;
		}
		amvTime = ArmV(time1);
		amvTime = AvgTrustedInterval(amvTime, time1);
		b << amvTime << endl;
		amvTime = ArmV(time2);
		amvTime = AvgTrustedInterval(amvTime, time2);
		g << amvTime << endl;
		amvTime = ArmV(time3);
		amvTime = AvgTrustedInterval(amvTime, time3);
		fR << amvTime << endl;
		amvTime = ArmV(time4);
		amvTime = AvgTrustedInterval(amvTime, time4);
		fp << amvTime << endl;
		amvTime = ArmV(time5);
		amvTime = AvgTrustedInterval(amvTime, time5);
		dp << amvTime << endl;
		amvTime = ArmV(time);
		amvTime = AvgTrustedInterval(amvTime, time);
		//cout << "Параллельно ";
		cout << amvTime << endl;
		string name = "c:\\test2\\par"  + to_string(count) + to_string(width) + "x" + to_string(height) + ".bmp";
		BMPWrite(rgb, header, bmiHeader, name.c_str());
		count++;
	}
	count = 2;
	while (count < 5)
	{
		tbb::global_control c(tbb::global_control::max_allowed_parallelism, count);
		for (int i = 0; i < 10; i++)
		{
			double t_start = omp_get_wtime();
			//Формирование массива яркости изображения
			double t = omp_get_wtime();
			brighti(height, width, rgb, I, Inew);
			double te = omp_get_wtime();
			time1[i] = (te - t) * 1000;
			//Вычисление первого и второго градиента яркости по горизонтали и вертикали
			t = omp_get_wtime();
			gradienti(height, width, Inew, DiffX, DiffY, DiffXY);
			te = omp_get_wtime();
			time2[i] = (te - t) * 1000;
			//Вычисление значения функции отклика угла R
			t = omp_get_wtime();
			funcRi(height, width, Inew, H, R, DiffX, DiffY, DiffXY, 20.0);
			te = omp_get_wtime();
			time3[i] = (te - t) * 1000;
			//Поиск особых точек для визуализации
			t = omp_get_wtime();
			findpointsi(height, width, H, HP, 20);
			te = omp_get_wtime();
			time4[i] = (te - t) * 1000;
			//Визуализация особых точек
			t = omp_get_wtime();
			drawpointsi(height, width, HP, rgb);
			te = omp_get_wtime();
			time5[i] = (te - t) * 1000;
			double t_end = omp_get_wtime();
			time[i] = (t_end - t_start) * 1000;
		}
		amvTime = ArmV(time1);
		amvTime = AvgTrustedInterval(amvTime, time1);
		b << amvTime << endl;
		amvTime = ArmV(time2);
		amvTime = AvgTrustedInterval(amvTime, time2);
		g << amvTime << endl;
		amvTime = ArmV(time3);
		amvTime = AvgTrustedInterval(amvTime, time3);
		fR << amvTime << endl;
		amvTime = ArmV(time4);
		amvTime = AvgTrustedInterval(amvTime, time4);
		fp << amvTime << endl;
		amvTime = ArmV(time5);
		amvTime = AvgTrustedInterval(amvTime, time5);
		dp << amvTime << endl;
		amvTime = ArmV(time);
		amvTime = AvgTrustedInterval(amvTime, time);
		//cout << "intel ";
		cout << amvTime << endl;
		string name = "c:\\test2\\intel" + to_string(count) + to_string(width) + "x" + to_string(height) + ".bmp";
		BMPWrite(rgb, header, bmiHeader, name.c_str());
		count++;
	}
	b.close(); g.close(); fR.close(); fp.close(); dp.close();
}
int main()
{
	ofstream b, g, fR, fp, dp;
	b.open("c:\\test2\\bright.txt", ios::trunc);
	g.open("c:\\test2\\g.txt", ios::trunc);
	fR.open("c:\\test2\\fr.txt", ios::trunc);
	fp.open("c:\\test2\\fp.txt", ios::trunc);
	dp.open("c:\\test2\\dp.txt", ios::trunc);
	b.close(); g.close(); fR.close(); fp.close(); dp.close();
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);
	RGBQUAD** rgb1;
	BITMAPFILEHEADER header1;
	BITMAPINFOHEADER bmiHeader1;
	BMPRead(rgb1, header1, bmiHeader1, "c:\\test2\\640x480.bmp");
	RGBQUAD** rgb2;
	BITMAPFILEHEADER header2;
	BITMAPINFOHEADER bmiHeader2;
	BMPRead(rgb2, header2, bmiHeader2, "c:\\test2\\1280x800.bmp");
	RGBQUAD** rgb3;
	BITMAPFILEHEADER header3;
	BITMAPINFOHEADER bmiHeader3;
	BMPRead(rgb3, header3, bmiHeader3, "c:\\test2\\1680x1050.bmp");
	RGBQUAD** rgb4;
	BITMAPFILEHEADER header4;
	BITMAPINFOHEADER bmiHeader4;
	BMPRead(rgb4, header4, bmiHeader4, "c:\\test2\\1920x1200.bmp");
	RGBQUAD** rgbs[] = { rgb1, rgb2, rgb3, rgb4 };
	BITMAPFILEHEADER headers[] = { header1, header2,  header3,  header4 };
	BITMAPINFOHEADER bmiHeaders[] = { bmiHeader1, bmiHeader2, bmiHeader3, bmiHeader4 };
	string name[4]{ "640x480", "1280x800", "1680x1050", "1920x1200" };
	for (int i = 0; i < 4; i++)
	{
		cout << name[i] << endl;
		formpictures(rgbs[i], bmiHeaders[i].biWidth, bmiHeaders[i].biHeight, headers[i], bmiHeaders[i]);
	}
	system("pause");
	return 0;
}
