#pragma once

//相邻特诊点结构
struct neighbor {
	int x;//x坐标
	int y;//y坐标
	int type;//1-端点，3-分岔点
	float Theta;//两点连线角度
	float Theta2Ridge; //两点脊线方向夹角
	float TheaThisNibor;//相邻特征点的脊线方向
	int distance;//两点距离
};

//特征点结构
struct minutiae {
	int x;
	int y;
	int type;
	float theta;//脊线方向
	neighbor * neibors;//相邻点特征信息
};

class FingerSolve
{
public:
	FingerSolve();
	~FingerSolve();

	//图像处理函数
	static int readBMP2buf(char * filename, unsigned char ** data, int &width, int&height, int&iDepth);
	static int MidFilter(unsigned char *ucImg, unsigned char *ucDstImg, int iWidth, int iHeight);
	static void Sort(unsigned char *data, int dsize);
	static int HistoNormalize(unsigned char*ucImg, unsigned char*ucNormImg, int iWidth, int iHeight);
	static 	int ImgDirection(unsigned char *ucImg, float *fDirc, int iWidth, int iHeight);
	static int DircLowPass(float *fDirc, float *fFitDirc, int iWidth, int iHeight);
	static int Frequency(unsigned char * ucImg, float * fDirection, float *fFrequency, int iWidth, int iHeight);
	static int GetMask(unsigned char *unImg, float *fDirection, float *fFrequency, unsigned char *ucMask, int iWidth, int iHeight);
	static int GaborEnhance(unsigned char *ucImg, float*fDirection, float *fFrequency, unsigned char *ucMask, unsigned char *ucImgEnhanced, int iWidth, int iHeight);
	static int BinaryImg(unsigned char * ucImage, unsigned char * ucBinImage, int iWidth, int iHeight, unsigned char uThreshold);
	static int BinaryToGray(unsigned char *ucBinImg, unsigned char *ucGrayImg, int iWidth, int iHeight);
	static int Thinning(unsigned char * ucBinedImg, unsigned char * ucThinnedImage,
		int iWidth, int iHeight, int iIterativeLimit);
	static int Extract(unsigned char *ucThinImg, unsigned char *ucMinuImg, int iWidth, int iHeight);
	static int DeEdgeMinu(minutiae * minutiaes, int count,unsigned char * ucImg,int iWidth,int iHeight);
	static int MinuFilter(byte * minuData, byte * thinData, minutiae *minutiaes, int &minuCount, int iWidth, int iHeight);


	static int SaveDataToImageFile(char *srcFile, char *dstFile, unsigned char *data);
	static int SaveDataToImageFile(char *srcFile, char *dstFile, float *data, int iWidth, int iHeight, int iDepth, float scale);


	static int WriteBMPImgFile(char *dstFileNmae, unsigned char **pusImgData);
	static wchar_t * ToWideChar(char * str);


};

