#include "stdafx.h"
#include "FingerSolve.h"
#include <conio.h>


FingerSolve::FingerSolve()
{
}


FingerSolve::~FingerSolve()
{
}

//写入bmp图片
int FingerSolve::WriteBMPImgFile(char * dstFileName, unsigned char ** pusImgData)
{
	FILE *fp;
	fopen_s(&fp, dstFileName, "r+b");
	if (!fp)
	{
		return -1;
	}
	int imgType, iWidth, iHeight;
	int iStartPos = 0;
	fseek(fp, 10L, SEEK_SET);
	fread((char *)(&iStartPos), 4, 1, fp);
	fseek(fp, 18L, SEEK_SET);
	fread((char*)(&iWidth), 4, 1, fp);
	fread((char*)(&iHeight), 4, 1, fp);
	unsigned short temp;
	fseek(fp, 28L, SEEK_SET);
	fread((char*)(&temp), 2, 1, fp);
	imgType = temp;
	if (imgType != 8)
	{
		return -2;
	}
	unsigned char*usImgData = *pusImgData;
	int iWidthInFile = 0;
	if (iWidth % 4 > 0)
	{
		iWidthInFile = iWidth - iWidth % 4 + 4;
	}
	else
	{
		iWidthInFile = iWidth;
	}
	for (int i = iHeight - 1; i >= 0; i--)
	{
		fseek(fp, iStartPos, SEEK_SET);
		fwrite((usImgData + i * iWidth), 1, iWidth, fp);
		iStartPos += iWidthInFile;

	}
	fclose(fp);
	return 0;
}
wchar_t * FingerSolve::ToWideChar(char * str)
{
	int num = MultiByteToWideChar(0, 0, str, -1, NULL, 0);
	wchar_t *wideStr = new wchar_t[num];
	MultiByteToWideChar(0, 0, str, -1, wideStr, num);
	return wideStr;
}

//从bmp文件读取位图数据
int FingerSolve::readBMP2buf(char * filename, unsigned char ** data, int &rWidth, int&rHeight, int&iDepth)
{
	CImage image;

	HRESULT r = image.Load(ToWideChar(filename));
	if (FAILED(r) || image.IsNull()) {
		return -1;
	}
	int width = image.GetWidth();
	int height = image.GetHeight();
	int bitDeep = image.GetBPP();//获取位图深度
	int pitch = image.GetPitch();//每行像素的长度

	if (bitDeep != 8)
	{
		return -1;
	}

	//	data = new unsigned char[pitch*height];//申请内存空间
	*data = (byte*)malloc((sizeof(byte)*width*height));
	memset(*data, 0, sizeof(byte)*width*height);

	byte *pData = NULL;
	if (image.GetPitch() < 0) //需要判断是自上向下还是自下向上的位图
		pData = (byte*)image.GetBits() + (image.GetPitch()*(image.GetHeight() - 1));
	else
		pData = (byte*)image.GetBits();

	byte gray = 0;

	byte *pRowData, *pRowImage, *pPixData, *pPixImage;

	//for (int y = 0; y < height; y++)
	//{
	//	pRowData = *data + width * y;
	//	pRowImage = pData + width * y;
	//	for (int x = 0; x < width; x++)
	//	{
	//		pPixData = pRowData + x;
	//		pPixImage = pRowImage + x;
	//		*pPixData = *pPixImage;
	//	}
	//}
	memcpy(*data, pData, sizeof(unsigned char)*width*height);

	rWidth = width;
	rHeight = height;

	image.Destroy();

	return 0;
}


//中值滤波
int  FingerSolve::MidFilter(unsigned char *ucImg, unsigned char *ucDstImg, int iWidth, int iHeight)
{
	memset(ucDstImg, 0, iWidth*iHeight);
	unsigned char *pUp, *pDown, *pImg;
	unsigned char x[9];
	for (int i = 1; i < iHeight - 1; i++)
	{
		pUp = ucImg + (i - 1)*iWidth;
		pImg = ucImg + i * iWidth;
		pDown = ucImg + (i + 1)*iWidth;


		for (int j = 1; j < iWidth - 1; j++)
		{
			pUp++;
			pImg++;
			pDown++;

			x[0] = *(pUp - 1);
			x[1] = *(pImg - 1);
			x[2] = *(pDown - 1);

			x[3] = *pUp;
			x[4] = *pImg;
			x[5] = *pDown;


			x[6] = *(pUp + 1);
			x[7] = *(pImg + 1);
			x[8] = *(pDown + 1);

			Sort(x, 9);
			*(ucDstImg + i * iWidth + j) = x[5];

		}
	}
	pDown = ucImg + iWidth;
	for (int j = 1; j < iWidth - 1; j++)
	{
		x[0] = *(ucImg + j - 1);
		x[1] = *(ucImg + j);
		x[2] = *(ucImg + j + 1);


		x[3] = *(pDown + j - 1);
		x[4] = *(pDown + j);
		x[5] = *(pDown + j + 1);


		Sort(x, 6);

		*(ucDstImg + j) = x[3];
		//*(ucDstImg + j) = x[0];

	}
	pUp = ucImg + iWidth * (iHeight - 2);
	pDown = ucImg + iWidth * (iHeight - 1);
	for (int j = 1; j < iWidth; j++)
	{
		x[0] = *(pDown + j - 1);
		x[1] = *(pDown + j);
		x[2] = *(pDown + j + 1);

		x[3] = *(pUp + j - 1);
		x[4] = *(pUp + j);
		x[5] = *(pUp + j + 1);

		Sort(x, 6);

		*(ucDstImg + iWidth * (iHeight - 1) + j) = x[3];
		//*(ucDstImg + iWidth * (iHeight - 1) + j) = x[0];

	}

	x[0] = *(ucImg);
	x[1] = *(ucImg + 1);
	x[2] = *(ucImg + iWidth);
	x[3] = *(ucImg + iWidth + 1);

	Sort(x, 4);


	*(ucDstImg) = x[2];

	x[0] = *(ucImg + iWidth - 1);
	x[1] = *(ucImg + iWidth - 2);
	x[2] = *(ucImg + 2 * iWidth - 1);
	x[3] = *(ucImg + 2 * iWidth - 2);


	Sort(x, 4);
	*(ucDstImg + (iHeight - 0)*iWidth - 1) = x[2];

	x[0] = *(ucImg + iHeight - 1);
	x[1] = *(ucImg + iHeight - 2);
	x[2] = *(ucImg + (iHeight - 1) * iWidth + 1);
	x[3] = *(ucImg + (iHeight - 2) * iWidth + 1);


	Sort(x, 4);
	*(ucDstImg + (iHeight - 0)*iWidth - 1) = x[2];



	return 0;


}


//排序
void FingerSolve::Sort(unsigned char * data, int dsize)
{
	unsigned char temp = 0;
	for (int i = 1; i < dsize; i++)
	{
		for (int j = dsize - 1; j > i; j--)
		{
			if (data[j] < data[i - 1])
			{
				temp = data[i];
				data[j] = data[j - 1];
				data[j - 1] = temp;
			}

		}

	}

}

//直方图均衡化
int FingerSolve::HistoNormalize(unsigned char * ucImg, unsigned char * ucNormImg, int iWidth, int iHeight)
{
	unsigned int Histogram[256];
	memset(Histogram, 0, 256 * sizeof(int));
	for (int i = 0; i < iHeight; i++)
	{
		for (int j = 0; j < iWidth; j++)
		{
			Histogram[ucImg[i*iWidth + j]]++;
		}

	}
	double dMean = 0;
	for (int i = 1; i < 255; i++)
	{
		dMean += i * Histogram[i];
	}
	dMean = int(dMean / (iWidth*iHeight));
	double dSigma = 0;
	for (int i = 0; i < 255; i++)
	{
		dSigma += Histogram[i] * (i - dMean)*(i - dMean);

	}
	dSigma /= (iWidth*iHeight);
	dSigma = sqrt(dSigma);
	double dMean0 = 128, dSigma0 = 128;
	double dCoeff = dSigma0 / dSigma;
	for (int i = 0; i < iHeight; i++)
	{
		for (int j = 0; j < iWidth; j++)
		{
			double dVal = ucImg[i*iWidth + j];
			dVal = dMean0 + dCoeff * (dVal - dMean0);
			if (dVal < 0)
			{

				dVal = 0;
			}
			else if (dVal > 255)
			{
				dVal = 255;

			}
			ucNormImg[i*iWidth + j] = (unsigned char)dVal;


		}

	}
	return 0;

}

//指纹脊线方向计算
int FingerSolve::ImgDirection(unsigned char * ucImg, float * fDirc, int iWidth, int iHeight)
{
	const int SEMISIZ = 7;
	int dx[SEMISIZ * 2 + 1][SEMISIZ * 2 + 1];
	int dy[SEMISIZ * 2 + 1][SEMISIZ * 2 + 1];
	float fx, fy;
	memset(fDirc, 0, iWidth*iHeight * sizeof(float));
	for (int y = SEMISIZ + 1; y < iHeight - SEMISIZ - 1; y++)
	{
		for (int x = SEMISIZ + 1; x < iWidth - SEMISIZ - 1; x++)
		{
			for (int j = 0; j < SEMISIZ * 2 + 1; j++)
			{
				for (int i = 0; i < SEMISIZ * 2 + 1; i++)
				{
					int index1 = (y + j - SEMISIZ)*iWidth + x + i - SEMISIZ;
					int index2 = (y + j - SEMISIZ)*iWidth + x + i - SEMISIZ - 1;
					int index3 = (y + j - SEMISIZ - 1)*iWidth + x + i - SEMISIZ;
					dx[i][j] = int(ucImg[index1] - ucImg[index2]);
					dy[i][j] = int(ucImg[index1] - ucImg[index3]);

					if (abs(dx[i][j]) > 50) {
						if (dx[i][j] > 0)
						{
							dx[i][j] = 255;
						}
						else
						{
							dx[i][j] = -255;
						}
					}
					else
					{
						dx[i][j] = 0;
					}

					if (abs(dy[i][j]) > 50) {
						if (dy[i][j] > 0)
						{
							dy[i][j] = 255;
						}
						else
						{
							dy[i][j] = -255;
						}
					}
					else
					{
						dy[i][j] = 0;
					}
				}
			}


			fx = 0.0;
			fy = 0.0;
			for (int j = 0; j < SEMISIZ * 2 + 1; j++)
			{
				for (int i = 0; i < SEMISIZ * 2 + 1; i++)
				{
					fx += 2 * dx[i][j] * dy[i][j];
					fy += (dx[i][j] * dx[i][j] - dy[i][j] * dy[i][j]);
				}

			}

			fDirc[y*iWidth + x] = atan2(fx, fy);

		}
	}
	return 0;

}

//低通滤波
int FingerSolve::DircLowPass(float * fDirc, float * fFitDirc, int iWidth, int iHeight)
{
	const int DIR_FILTER_SIZE = 2;
	int blocksize = 2 * DIR_FILTER_SIZE + 1;
	int imgsize = iWidth * iHeight;

	float *filter = new float[blocksize*blocksize];
	float *phix = new float[imgsize];
	float *phiy = new float[imgsize];
	float *phi2x = new float[imgsize];
	float *phi2y = new float[imgsize];

	memset(fFitDirc, 0, sizeof(float)*iWidth*iHeight);

	float  tempSum = 0.0;
	for (int y = 0; y < blocksize; y++)
	{

		for (int x = 0; x < blocksize; x++)
		{

			filter[y*blocksize + x] = (float)(blocksize - (abs(DIR_FILTER_SIZE - x) + abs(DIR_FILTER_SIZE - y)));
			tempSum += filter[y*blocksize + x];
		}
	}

	for (int y = 0; y < blocksize; y++)
	{
		for (int x = 0; x < blocksize; x++)
		{
			filter[y*blocksize + x] /= tempSum;
		}
	}

	for (int y = 0; y < iHeight; y++)
	{
		for (int x = 0; x < iWidth; x++)
		{
			phix[y*iWidth + x] = cos(fDirc[y*iWidth + x]);
			phiy[y*iWidth + x] = sin(fDirc[y*iWidth + x]);

		}
	}
	memset(phi2x, 0, sizeof(float)*imgsize);
	memset(phi2y, 0, sizeof(float)*imgsize);
	float nx, ny;
	int val;
	for (int y = 0; y < iHeight - blocksize; y++)
	{
		for (int x = 0; x < iWidth - blocksize; x++)
		{
			nx = 0.0;
			ny = 0.0;
			for (int j = 0; j < blocksize; j++)
			{
				for (int i = 0; i < blocksize; i++)
				{


					val = (x + i) + (j + y)*iWidth;
					nx += filter[j*blocksize + i] * phix[val];
					ny += filter[j*blocksize + i] * phiy[val];
				}
			}
			val = x + y * iWidth;
			phi2x[val] = nx;
			phi2y[val] = ny;


		}

	}
	for (int y = 0; y < iHeight - blocksize; y++)
	{
		for (int x = 0; x < iWidth - blocksize; x++)
		{

			val = x + y * iWidth;
			fFitDirc[val] = atan2(phi2y[val], phi2x[val])*0.5;

		}
	}
	delete[] phi2y;
	delete[] phi2x;
	delete[] phiy;
	delete[] phix;

	return 0;
}

//频率计算
int FingerSolve::Frequency(unsigned char * ucImg, float * fDirection, float * fFrequency, int iWidth, int iHeight)
{
	const int SIZE_L = 32;
	const int SIZE_W = 16;
	const int SIZE_L2 = 16;
	const int SIZE_W2 = 8;

	int peak_pos[SIZE_L];
	int peak_cnt;
	float peak_freq;
	float Xsig[SIZE_L];

	float dir = 0.0;
	float cosdir = 0.0;
	float sindir = 0.0;
	float maxPeak, minPeak;

	float *frequency1 = new float[iWidth*iHeight];
	memset(fFrequency, 0, sizeof(float)*iWidth*iHeight);
	memset(frequency1, 0, sizeof(float)*iWidth*iHeight);

	int x, y;
	int d, k;
	int u, v;

	for (y = SIZE_L2; y < iHeight - SIZE_L2; y++)
	{
		for (x = SIZE_L2; x < iWidth - SIZE_L2; x++)
		{
			dir = fDirection[(y + SIZE_W2)*iWidth + (x + SIZE_W2)];
			cosdir = -sin(dir);
			sindir = cos(dir);

			for (k = 0; k < SIZE_L; k++)
			{
				Xsig[k] = 0.0;
				for (d = 0; d < SIZE_W; d++)
				{
					u = (int)(x + (d - SIZE_W2)*cosdir + (k - SIZE_L2)*sindir);
					v = (int)(y + (d - SIZE_W2)*sindir - (k - SIZE_L2)*cosdir);
					if (u < 0)
					{
						u = 0;
					}
					else if (u > iWidth - 1)
					{
						u = iWidth - 1;
					}
					if (v < 0)
					{
						v = 0;
					}
					else if (v > iHeight - 1)
					{
						v = iHeight - 1;
					}

					Xsig[k] += ucImg[u + v * iWidth];
				}
				Xsig[k] /= SIZE_W;
			}
			maxPeak = minPeak = Xsig[0];
			for (k = 0; k < SIZE_L; k++)
			{
				if (minPeak > Xsig[k])
				{
					minPeak = Xsig[k];
				}
				if (maxPeak < Xsig[k])
				{
					maxPeak = Xsig[k];
				}
			}

			peak_cnt = 0;
			if ((maxPeak - minPeak) > 64)
			{
				for (k = 0; k < SIZE_L; k++)
				{
					if ((Xsig[k - 1] < Xsig[k]) && (Xsig[k] >= Xsig[k + 1]))
					{
						peak_pos[peak_cnt++] = k;
					}
				}
			}

			peak_freq = 0.0;
			if (peak_cnt >= 2)
			{
				for (k = 0; k < peak_cnt - 1; k++)
				{
					peak_freq += (peak_pos[k + 1] - peak_pos[k]);
				}
				peak_freq /= peak_cnt - 1;
			}
			if (peak_freq<3.0 || peak_freq>25.0)
			{
				frequency1[x + y * iWidth] = 0.0;
			}
			else
			{
				frequency1[x + y * iWidth] = 1.0 / peak_freq;
			}
		}
	}

	for (y = SIZE_L2; y < iHeight - SIZE_L2; y++)
	{
		for (x = SIZE_L2; x < iWidth - SIZE_L2; x++)
		{
			k = x + y * iWidth;
			peak_freq = 0.0;
			for (v = -2; v <= 2; v++)
			{
				for (u = -2; u <= 2; u++)
				{
					peak_freq += frequency1[(x + u) + (y + v)*iWidth];
				}
			}
			fFrequency[k] = peak_freq / 25;
		}
	}

	delete[] frequency1;

	return 0;
}

//掩码计算
int FingerSolve::GetMask(unsigned char * unImg, float * fDirection, float * fFrequency, unsigned char * ucMask, int iWidth, int iHeight)
{
	float freqMin = 1.0 / 25.0;
	float freqMax = 1.0 / 3.0;
	int x, y, k;
	int pos, posout;

	memset(ucMask, 0, iWidth*iHeight);

	for (y = 0; y < iHeight; y++)
	{
		for (x = 0; x < iWidth; x++)
		{
			pos = x + y * iWidth;
			posout = x + y * iWidth;
			ucMask[posout] = 0;
			if (fFrequency[pos] >= freqMin && fFrequency[pos] <= freqMax)
			{
				ucMask[posout] = 255;

			}
		}
	}

	for (k = 0; k < 4; k++)
	{
		for (y = 1; y < iHeight - 1; y++)
		{
			for (x = 1; x < iWidth - 1; x++)
			{
				if (ucMask[x + y * iWidth] == 0xFF)
				{
					ucMask[x - 1 + y * iWidth] |= 0x80;
					ucMask[x + 1 + y * iWidth] |= 0x80;
					ucMask[x + (y - 1)* iWidth] |= 0x80;
					ucMask[x + (y + 1)* iWidth] |= 0x80;
				}
			}
		}
		for (y = 1; y < iHeight - 1; y++)
		{
			for (x = 1; x < iWidth - 1; x++)
			{
				if (ucMask[x + y * iWidth])
				{
					ucMask[x + y * iWidth] = 0xFF;
				}
			}
		}


	}

	for (k = 0; k < 12; k++)
	{
		for (y = 1; y < iHeight - 1; y++)
		{
			for (x = 1; x < iWidth - 1; x++)
			{
				if (ucMask[x + y * iWidth] == 0x0)
				{
					ucMask[x - 1 + y * iWidth] &= 0x80;
					ucMask[x + 1 + y * iWidth] &= 0x80;
					ucMask[x + (y - 1) * iWidth] &= 0x80;
					ucMask[x + (y + 1) * iWidth] &= 0x80;
				}
			}
		}
		for (y = 1; y < iHeight - 1; y++)
		{
			for (x = 1; x < iWidth - 1; x++)
			{
				if (ucMask[x + y * iWidth] != 0xFF)
				{
					ucMask[x + y * iWidth] = 0x0;
				}
			}
		}
	}
	return 0;
}

//gabor增强
int FingerSolve::GaborEnhance(unsigned char * ucImg, float*fDirection, float * fFrequency, unsigned char * ucMask, unsigned char * ucImgEnhanced, int iWidth, int iHeight)
{
	const float PI = 3.141592654;
	int i, j, u, v;
	int wg2 = 5;
	float sum, f, g;
	float x2, y2;
	float dx2 = 1.0*(4.0*4.0);
	float dy2 = 1.0*(4.0*4.0);

	memset(ucImgEnhanced, 0, iWidth*iHeight);


	for (j = wg2; j < iHeight - wg2; j++)
	{
		for (i = wg2; i < iWidth - wg2; i++)
		{

			if (ucMask[i + j * iWidth] == 0)
			{
				continue;
			}
			g = fDirection[i + j * iWidth];
			f = fFrequency[i + j * iWidth];
			g += PI / 2;
			sum = 0.0;
			for (v = -wg2; v <= wg2; v++)
			{
				for (u = -wg2; u <= wg2; u++)
				{
					x2 = -u * sin(g) + v * cos(g);
					y2 = u * cos(g) + v * sin(g);
					sum += exp(-0.5*(x2*x2*dx2 + y2 * y2*dy2))*cos(2 * PI*x2*f)*ucImg[(i - u) + (j - v)*iWidth];
				}
			}

			if (sum > 255.0)
			{
				sum = 255.0;
			}
			if (sum < 0.0)
			{
				sum = 0.0;
			}
			ucImgEnhanced[i + j * iWidth] = (unsigned char)sum;

		}
	}
	return 0;

}

//计算二值化图像
int FingerSolve::BinaryImg(unsigned char * ucImage, unsigned char * ucBinImage, int iWidth, int iHeight, unsigned char uThreshold)
{
	unsigned char *pStart = ucImage, *pEnd = ucImage + iWidth * iHeight;
	unsigned char *pDest = ucBinImage;
	while (pStart < pEnd)
	{
		*pDest = *pStart > uThreshold ? 1 : 0;
		pStart++;
		pDest++;
	}
	return 0;
}

//转化为256灰度图像
int FingerSolve::BinaryToGray(unsigned char * ucBinImg, unsigned char * ucGrayImg, int iWidth, int iHeight)
{
	unsigned char *pStart = ucBinImg, *pEnd = ucBinImg + iWidth * iHeight;
	unsigned char *pDest = ucGrayImg;

	while (pStart < pEnd)
	{
		*pDest = (*pStart) > 0 ? 255 : 0;
		pStart++;
		pDest++;
	}

	return 0;
}

//指纹细化
int FingerSolve::Thinning(unsigned char * ucBinedImg, unsigned char * ucThinnedImage, int iWidth, int iHeight, int iIterativeLimit)
{
	unsigned char x1, x2, x3, x4, x5, x6, x7, x8, xp;
	unsigned char g1, g2, g3, g4;
	unsigned char b1, b2, b3, b4;
	unsigned char np1, np2, npm;
	unsigned char *pUp, *pDown, *pImg;
	int iDeletedPoints = 0;

	memcpy(ucThinnedImage, ucBinedImg, iWidth*iHeight);

	for (int it = 0; it < iIterativeLimit; it++)
	{
		iDeletedPoints = 0;

		for (int i = 1; i < iHeight - 1; i++)
		{
			pUp = ucBinedImg + (i - 1)*iWidth;
			pImg = ucBinedImg + i * iWidth;
			pDown = ucBinedImg + (i + 1)*iWidth;

			for (int j = 1; j < iWidth - 1; j++)
			{
				pUp++;
				pImg++;
				pDown++;

				if (!*pImg)
				{
					continue;
				}

				x6 = *(pUp - 1);
				x5 = *(pImg - 1);
				x4 = *(pDown - 1);
				x7 = *pUp;
				xp = *pImg;
				x3 = *pDown;

				x8 = *(pUp + 1);
				x1 = *(pImg + 1);
				x2 = *(pDown + 1);

				b1 = !x1 && (x2 == 1 || x3 == 1);
				b2 = !x3 && (x4 == 1 || x5 == 1);
				b3 = !x5 && (x6 == 1 || x7 == 1);
				b4 = !x7 && (x8 == 1 || x1 == 1);

				g1 = (b1 + b2 + b3 + b4) == 1;

				np1 = x1 || x2;
				np1 += x3 || x4;
				np1 += x5 || x6;
				np1 += x7 || x8;
				np2 = x2 || x3;
				np2 += x4 || x5;
				np2 += x6 || x7;
				np2 += x8 || x1;

				npm = np1 > np2 ? np2 : np1;
				g2 = npm >= 2 && npm <= 3;
				g3 = (x1 && (x2 || x3 || !x8)) == 0;
				g4 = (x5 && (x6 || x7 || !x4)) == 0;

				if (g1&&g2&&g3)
				{
					ucThinnedImage[iWidth*i + j] = 0;
					++iDeletedPoints;
				}
			}
		}

		memcpy(ucBinedImg, ucThinnedImage, iWidth*iHeight);

		for (int i = 1; i < iHeight - 1; i++)
		{
			pUp = ucBinedImg + (i - 1)*iWidth;
			pImg = ucBinedImg + i * iWidth;
			pDown = ucBinedImg + (i + 1)*iWidth;

			for (int j = 1; j < iWidth - 1; j++)
			{
				pUp++;
				pImg++;
				pDown++;
				if (!*pImg)
				{
					continue;
				}

				x6 = *(pUp - 1);
				x5 = *(pImg - 1);
				x4 = *(pDown - 1);

				x7 = *pUp;
				xp = *pImg;
				x3 = *pDown;

				x8 = *(pUp + 1);
				x1 = *(pImg + 1);
				x2 = *(pDown + 1);

				b1 = !x1 && (x2 == 1 || x3 == 1);
				b2 = !x3 && (x4 == 1 || x5 == 1);
				b3 = !x5 && (x6 == 1 || x7 == 1);
				b4 = !x7 && (x8 == 1 || x1 == 1);

				g1 = (b1 + b2 + b3 + b4) == 1;

				np1 = x1 || x2;
				np1 += x3 || x4;
				np1 += x5 || x6;
				np1 += x7 || x8;
				np2 = x2 || x3;
				np2 += x4 || x5;
				np2 += x6 || x7;
				np2 += x8 || x1;

				npm = np1 > np2 ? np2 : np1;
				g2 = npm >= 2 && npm <= 3;
				g3 = (x1 && (x2 || x3 || !x8)) == 0;
				g4 = (x5 && (x6 || x7 || !x4)) == 0;

				if (g1&&g2&&g4)
				{
					ucThinnedImage[iWidth*i + j] = 0;
					++iDeletedPoints;
				}
			}
		}
		memcpy(ucBinedImg, ucThinnedImage, iWidth*iHeight);

		if (iDeletedPoints == 0)
		{
			break;
		}
	}

	for (int i = 0; i < iHeight; i++)
	{
		for (int j = 0; j < iWidth; j++)
		{
			if (i < 16)
			{
				ucThinnedImage[i*iWidth + j] = 0;
			}
			else if (i >= iHeight - 16)
			{
				ucThinnedImage[i*iWidth + j] = 0;
			}
			else if (j < 16)
			{
				ucThinnedImage[i*iWidth + j] = 0;
			}
			else if (j >= iWidth - 16)
			{
				ucThinnedImage[i*iWidth + j] = 0;
			}
		}
	}
	return 0;
}

//指纹特征提取
int FingerSolve::Extract(unsigned char * ucThinImg, unsigned char * ucMinuImg, int iWidth, int iHeight)
{
	unsigned char *pDest = ucMinuImg;
	unsigned char *pUp, *pImg, *pDown;
	unsigned char x1, x2, x3, x4, x5, x6, x7, x8;
	unsigned char nc;
	int iMinuCount = 0;

	memset(pDest, 0, iWidth*iHeight);

	for (int i = 1; i < iHeight - 1; i++)
	{
		pUp = ucThinImg + (i - 1)*iWidth;
		pImg = ucThinImg + i * iWidth;
		pDown = ucThinImg + (i + 1)*iWidth;

		for (int j = 1; j < iWidth - 1; j++)
		{
			pUp++;
			pImg++;
			pDown++;

			if (!*pImg)
			{
				continue;
			}

			x6 = *(pUp - 1);
			x5 = *(pImg - 1);
			x4 = *(pDown - 1);

			x7 = *pUp;
			x3 = *pDown;

			x8 = *(pUp + 1);
			x1 = *(pImg + 1);
			x2 = *(pDown + 1);

			nc = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8;

			if (nc == 1)
			{
				pDest[i*iWidth + j] = 1;
				++iMinuCount;
			}
			else if (nc == 3)
			{
				pDest[i*iWidth + j] = 3;
				++iMinuCount;
			}
		}

	}
	return iMinuCount;

}

//特征点过滤
int FingerSolve::DeEdgeMinu(minutiae * minutiaes, int count, unsigned char * ucImg, int iWidth, int iHeight)
{
	int minuCount = count;
	int x, y, type;
	int center_x, center_y;
	bool del;
	double k;//斜率

	//初始化
	int *pFlag = new int[minuCount];
	k = iHeight * 1.0 / iWidth;
	//寻找指纹图像中心点
	center_x = 0, center_y = 0;
	for (int i = 0; i < minuCount; ++i) {
		center_x += minutiaes[i].x;
		center_y += minutiaes[i].y;
	}
	center_x /= minuCount;
	center_y /= minuCount;
	//舍弃指纹矩阵边角之外特征点
	int min_x = 0, max_x = 0, min_y = 0, max_y = 0;
	for (int i = 0; i < minuCount; ++i) {
		minutiaes[i].x > max_x ? max_x = minutiaes[i].x : 1;
		minutiaes[i].x < min_x ? min_x = minutiaes[i].x : 1;
		minutiaes[i].y < min_y ? min_y = minutiaes[i].y : 1;
		minutiaes[i].y > max_y ? max_y = minutiaes[i].y : 1;
	}

	double threshold = 0.25;
	//计算舍弃范围（离中心点在各方向上10%外的点）
	min_x = min_x + threshold * abs((center_x - min_x));
	max_x = max_x - threshold * abs((center_x - max_x));
	min_y = min_y + threshold * abs((center_y - min_y));
	max_y = max_y - threshold * abs((center_y - max_y));

	memset(pFlag, 0, sizeof(int)*minuCount);
	for (int i = 0; i < minuCount; ++i)
	{
		//当前特征点信息
		y = minutiaes[i].y;
		x = minutiaes[i].x;

		del = true;

		//if (y < -k * x + iHeight) {
		//	if (y < k*x)
		//	{
		//		//上
		//		while (--y>=0)
		//		{
		//			if (ucImg[x + y * iWidth] > 0) {
		//				del != del;
		//				break;
		//			}
		//		}
		//	}
		//	else
		//	{
		//		//左
		//		while (--x>=0)
		//		{
		//			if (ucImg[x + y*iWidth] > 0) {
		//				del != del;
		//				break;
		//			}
		//		}
		//	}
		//}
		//else
		//{
		//	if (y < k*x) {
		//		//右
		//		while (++x<iWidth)
		//		{
		//			if (ucImg[x + y * iWidth] > 0) {
		//				del != del;
		//				break;
		//			}
		//		}
		//	}
		//	else
		//	{
		//		//下
		//		while (++y<iHeight)
		//		{
		//			if (ucImg[x + y * iWidth] > 0) {
		//				del != del;
		//				break;
		//			}
		//		}
		//	}
		//}
		if (x > min_x&&x<max_x&&y>min_x&&y < max_y)
		{
			del = false;
		}

		//del==true -> 边缘点 -> 删除
		if (del)pFlag[i] = 1;//标记数组中设置为1（删除）

	}

	int newCount = 0;
	//特征重排
	for (int i = 0; i != minuCount; ++i)
	{
		if (pFlag[i] == 0)
		{
			memcpy(&minutiaes[newCount++], &minutiaes[i], sizeof(minutiae));
		}
	}

	delete[] pFlag;
	pFlag = NULL;

	return newCount;
}

int FingerSolve::MinuFilter(byte * minuData, byte * thinData, minutiae * minutiaes, int & minuCount, int iWidth, int iHeight)
{
	float *dirs = new float[iWidth*iHeight];

	memset(dirs, 0, iWidth*iHeight * sizeof(float));

	//计算脊线方向
	ImgDirection(thinData, dirs, iWidth, iHeight);

	//提取特征点
	byte * pImg;

	int tmp = 0;
	for (size_t i = 0; i < iHeight - 1; i++)
	{
		pImg = minuData + i * iWidth;
		for (size_t j = 0; j < iWidth - 1; j++)
		{
			++pImg;
			if (*pImg > 0)
			{
				minutiaes[tmp].x = j;
				minutiaes[tmp].y = i;
				minutiaes[tmp].theta = dirs[i*iWidth + j];
				minutiaes[tmp].type = *pImg;
				++tmp;
			}
		}
	}

	delete[]dirs;

	//去除边缘特征点
	minuCount = DeEdgeMinu(minutiaes, minuCount, thinData, iWidth, iHeight);

	//去除毛刺
	int *pFlag = new int[minuCount];
	memset(pFlag, 0, sizeof(int)*minuCount);

	int x1, x2, y1, y2, type1, type2;
	for (size_t i = 0; i < minuCount; i++)
	{
		x1 = minutiaes[i].x;
		y1 = minutiaes[i].y;
		type1 = minutiaes[i].type;
		for (size_t j = i + 1; j < minuCount; j++)
		{
			if (pFlag[j] == 1)continue;//已经删除的特征点,直接跳过

			x2 = minutiaes[j].x;
			y2 = minutiaes[j].y;
			type2 = minutiaes[j].type;

			int d = (int)sqrt(float((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2)));

			//临近点才有关系
			if (d <= 6)
			{
				if (type1 == type2)//两个点相同类型
				{
					if (type1 == 1) //间断
					{
						pFlag[i] = 1;
						pFlag[j] = 1;
					}
					else //小岛
					{
						pFlag[j] = 1;
					}
				}
				else if (type1 == 1) //点1是端点
				{
					pFlag[i] = 1;
				}
				else//点2是端点
				{
					pFlag[j] = 1;
				}
			}


		}
	}

	int newCount = 0;
	for (int i = 0; i != minuCount; ++i)
	{
		if (pFlag[i] == 0)
		{
			memcpy(&minutiaes[newCount++], &minutiaes[i], sizeof(minutiae));
		}
	}

	minuCount = newCount;

	delete[] pFlag;

	return 0;
}

//返回两个点的线段倾斜度
float FingerSolve::AngleOfPoints(int x1, int y1, int x2, int y2)
{
	const float PI = 3.1415925535;
	float deltaY, deltaX, theta = 0.0f;
	deltaX = x2 - x1;
	deltaY = y2 - y1;
	if (deltaX > 0 && deltaY < 0) {
		theta = atan2(deltaY, -1 * deltaX);
	}
	else if (deltaY < 0 && deltaX < 0)
	{
		theta = PI - atan2(-1 * deltaY, -1 * deltaX);
	}
	else if (deltaY < 0 && deltaX > 0)
	{
		theta = atan2(-1 * deltaY, deltaX);
	}
	else if (deltaY < 0 && deltaX < 0)
	{
		theta = PI - atan2(-1 * deltaY, -1 * deltaX);
	}
	else if (deltaX == 0)
	{
		theta = PI / 2;
	}
	else
	{
		theta = 0;
	}

	return theta;
}

//建立特征点之间的相邻关系
int FingerSolve::BuildNeighborsShip(minutiae * minus, int minuCount)
{
	const int MAX_NEIBORS = 10;//每个特征点最多的相邻点
	int x1, x2, y1, y2;
	int * pFlag = new int[minuCount];
	for (size_t i = 0; i < minuCount; i++)
	{
		x1 = minus[i].x;
		y1 = minus[i].y;

		//初始化当前特征点相邻特征点标记
		memset(pFlag, 0, sizeof(int)*minuCount);//全0（不相邻）
		pFlag[i] = 1;//将自身标记为相邻

		//初始化当前特征点相邻数组结构
		minus[i].neibors = new neighbor[MAX_NEIBORS];
		memset(minus[i].neibors, 0, sizeof(neighbor)*MAX_NEIBORS);

		for (size_t neiborNo = 0; neiborNo < MAX_NEIBORS; neiborNo++)
		{
			//查找距离最小的特征点
			unsigned int minDist = 0xffff;//最小距离特征点距离
			int minNo = 0;//最小距离特征点下标

			for (size_t j = 0; j < minuCount; j++)
			{
				if (pFlag[i] == 0)
				{
					x2 = minus[j].x;
					y2 = minus[j].y;
					int dtmp = (int)sqrt(float((y1 - y2)*(y1 - y2) + (x1 - x2)*(x1 - x2)));
					if (dtmp < minDist)
					{
						minDist = dtmp;
						minNo = j;
					}
				}
			}

			//保存结果
			pFlag[minNo] = 1;//设置相邻
			minus[i].neibors[neiborNo].x = minus[minNo].x;
			minus[i].neibors[neiborNo].y = minus[minNo].y;
			minus[i].neibors[neiborNo].type = minus[minNo].type;
			minus[i].neibors[neiborNo].Theta = AngleOfPoints(minus[minNo].x, minus[minNo].y, x1, y1);
			minus[i].neibors[neiborNo].Theta2Ridge = minus[minNo].theta - minus[minNo].theta;
			minus[i].neibors[neiborNo].TheaThisNibor = minus[minNo].theta;//相邻特征点的脊线方向
			minus[i].neibors[neiborNo].distance = minDist;
		}
	}

	delete[] pFlag;
	return 0;
}

//计算两个特征点的相似度
float FingerSolve::MinuSimilarity(int iWidth, int iHeight, minutiae * minutiae1, int count1, minutiae * minutiae2, int count2)
{
	const int MAX_SIMILAR_PAIR = 100;
	const int MAX_NEIGHBOR_EACH = 10;

	BuildNeighborsShip(minutiae1, count1);
	BuildNeighborsShip(minutiae2, count2);

	int similarPair[MAX_SIMILAR_PAIR][2];
	memset(similarPair, 0, 100 * 2 * sizeof(int));

	minutiae* baseminutiae;
	minutiae* refminutiae;
	int baseAccount, refAccount;

	if (count1 < count2)
	{
		baseminutiae = minutiae1;
		baseAccount = count1;
		refminutiae = minutiae2;
		refAccount = count2;
	}
	else
	{
		baseminutiae = minutiae2;
		baseAccount = count2;
		refminutiae = minutiae1;
		refAccount = count1;
	}

	neighbor *baseNeighbors = NULL;
	neighbor *refNeighbors = NULL;
	int similarminutiae = 0;
	float baseTheta, refTheta;

	for (int i = 0; i < baseAccount; i++)
	{
		baseNeighbors = baseminutiae[i].neibors;
		baseTheta = baseminutiae[i].theta;
		int refSimilarNo = 0;
		int maxSimilarNeibors = 0;
		for (int j = 0; j < refAccount; j++)
		{
			if (refminutiae[j].type != baseminutiae[i].type)
			{
				continue;
			}
			refNeighbors = refminutiae[j].neibors;
			refTheta = refminutiae[j].theta;


			int thisSimilarNeigbors = 0;

			for (int m = 0; m < MAX_NEIGHBOR_EACH; m++)
			{
				for (int n = 0; n < MAX_NEIGHBOR_EACH; n++)
				{
					if (baseNeighbors[m].type != refNeighbors[n].type)
					{
						continue;

					}
					int dist = abs(int(baseNeighbors[m].distance - refNeighbors[n].distance));
					float theta1 = fabs(float((baseNeighbors[m].Theta - baseTheta) - (refNeighbors[n].Theta - refTheta)));
					float theta2 = fabs(float(baseNeighbors[m].Theta2Ridge - refNeighbors[n].Theta2Ridge));
					float theta3 = fabs(float((baseNeighbors[m].Theta - baseNeighbors[m].TheaThisNibor) - (refNeighbors[n].Theta - refNeighbors[n].TheaThisNibor)));

					if (dist < 4 && theta1 < 0.15f&&theta2 < 0.15f&&theta3 < 0.15f);
					{
						++thisSimilarNeigbors;
						break;
					}
				}
			}
			if ((thisSimilarNeigbors >= MAX_NEIGHBOR_EACH * 3 / 10) && (similarminutiae < MAX_SIMILAR_PAIR))
			{
				similarPair[similarminutiae][0] = i;
				similarPair[similarminutiae][1] = refSimilarNo;
				++similarminutiae;
			}
		}
	}
	float similarity = similarminutiae / 8.0f;

	similarity = similarminutiae < 20 ? 0.0f : similarity;
	similarity = similarminutiae > 40 ? 1.0f : similarity;
	return similarity;
}

//保存图像
int FingerSolve::SaveDataToImageFile(char * srcFile, char * dstFile, unsigned char * data)
{
	CopyFile(ToWideChar(srcFile), ToWideChar(dstFile), false);
	WriteBMPImgFile(dstFile, &data);
	return 0;
}

//保存图像
int FingerSolve::SaveDataToImageFile(char * srcFile, char * dstFile, float * data, int iWidth, int iHeight, int iDepth, float scale)
{
	CopyFile(ToWideChar(srcFile), ToWideChar(dstFile), false);
	unsigned char *tmpData = new unsigned char[iWidth*iHeight];
	for (int i = 0; i<int(iWidth*iHeight); i++)
	{
		tmpData[i] = unsigned char((scale*data[i]));
	}
	WriteBMPImgFile(dstFile, &tmpData);

	delete[] tmpData;

	return 0;
}

