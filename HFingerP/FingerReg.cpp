// FingerReg.cpp: 实现文件
//

#include "stdafx.h"
#include "HFingerP.h"
#include "FingerReg.h"
#include "afxdialogex.h"

#define PATHCAT(a,b) a##b  
#define g(a) #a
#define h(a) g(a)

#ifdef _DEBUG
#define DEBUGDIR "..\\Debug\\"
#define DEBUGFINGERS "..\\..\\fingers_bmp\\"
#else
#define DEBUGDIR ".\\"
#define DEBUGFINGERS ".\\"
#endif



int ShowImageInCtrl1(CStatic &picCtrl, char*filename);

BOOL isRefresh = false;
BOOL isShowIng = false;

// FingerReg 对话框

IMPLEMENT_DYNAMIC(FingerReg, CDialog)

FingerReg::FingerReg(CWnd* pParent /*=nullptr*/)
	: CDialog(IDD_FINGERREG_DIALOG, pParent)
{
	flag = 1;
	mdb = new db(PATHCAT(DEBUGDIR,"miu.db"));
}

FingerReg::~FingerReg()
{
	
}


//初始化采集器
void FingerReg::initEngine()
{
	chgEgeState(3);
	long result = m_zkfpEng.InitEngine();//启动指纹采集器 result:0 连接成功，2 未连接指纹采集器
	if (result != 0) {
		chgEgeState(0);
		MessageBox(_T("指纹采集器未连接或故障！"), _T("提示"));
		return;
	}
	chgEgeState(1);
	/*(CButton*)this->GetDlgItem(IDC_BUTTON_INITENGENINE)->EnableWindow(false);
	(CButton*)this->GetDlgItem(IDC_BUTTON_STOPENGINE)->EnableWindow(true);
	(CButton*)this->GetDlgItem(IDC_BUTTON_REG)->EnableWindow(true);*/
}
void FingerReg::stopEngine()
{
	chgEgeState(4);
	m_zkfpEng.EndEngine();//断开指纹采集器
	chgEgeState(0);
	//(CButton*)this->GetDlgItem(IDC_BUTTON_REG)->EnableWindow(false);
	//(CButton*)this->GetDlgItem(IDC_BUTTON_STOPENGINE)->EnableWindow(false);
	//(CButton*)this->GetDlgItem(IDC_BUTTON_INITENGENINE)->EnableWindow(true);
	//MessageBox(_T("停止指纹采集器"), _T("提示"));
}

void FingerReg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_ZKFPENGX1, m_zkfpEng); // 绑定指纹采集器控件和变量
	DDX_Control(pDX, IDC_STATIC_PIC1, m_picCtrl); // 图像显示控件控件和变量
	//DDX_Control(pDX, IDC_BUTTON_INITENGENINE, m_btnInitEngine); // 启动指纹采集器和变量
	//DDX_Control(pDX, IDC_BUTTON_STOPENGINE, m_btnStopEngine); // 停止指纹采集器和变量
	DDX_Control(pDX, IDC_STATIC_TEXT_ENGINESTATE, m_textEngState); // 指纹采集器状态文字
	DDX_Control(pDX, IDC_STATIC_TEXT_ENGINESTATE, m_textEngState); // 文本框

}

BEGIN_EVENTSINK_MAP(FingerReg, CDialog)
	ON_EVENT(FingerReg, IDC_ZKFPENGX1, 8, OnImageReceivedZkfpengx, VTS_PBOOL)
END_EVENTSINK_MAP()

BEGIN_MESSAGE_MAP(FingerReg, CDialog)
	ON_BN_CLICKED(IDCANCEL, &FingerReg::OnBnClickedCancel)
	ON_BN_CLICKED(IDC_BUTTON_test, &FingerReg::OnBnClickedButtontest)
	/*添加RadioButton的消息处理*/
	ON_BN_CLICKED(IDC_RADIO_SOURCE, &FingerReg::OnBnClickedRadioShowImage)
	ON_BN_CLICKED(IDC_RADIO2, &FingerReg::OnBnClickedRadioShowImage)
	ON_BN_CLICKED(IDC_RADIO3, &FingerReg::OnBnClickedRadioShowImage)
	ON_BN_CLICKED(IDC_RADIO4, &FingerReg::OnBnClickedRadioShowImage)
	ON_BN_CLICKED(IDC_RADIO5, &FingerReg::OnBnClickedRadioShowImage)
	ON_BN_CLICKED(IDC_RADIO6, &FingerReg::OnBnClickedRadioShowImage)
	ON_BN_CLICKED(IDC_RADIO7, &FingerReg::OnBnClickedRadioShowImage)
	ON_BN_CLICKED(IDC_RADIO8, &FingerReg::OnBnClickedRadioShowImage)
	ON_BN_CLICKED(IDC_RADIO9, &FingerReg::OnBnClickedRadioShowImage)
	ON_BN_CLICKED(IDC_RADIO10, &FingerReg::OnBnClickedRadioShowImage)
	ON_BN_CLICKED(IDC_RADIO11, &FingerReg::OnBnClickedRadioShowImage)
	/*添加RadioButton的消息处理*/

	ON_BN_CLICKED(IDC_BUTTON_test_valid, &FingerReg::OnBnClickedButtontestvalid)
END_MESSAGE_MAP()


// FingerReg 消息处理程序

//退出对话框
void FingerReg::OnBnClickedCancel()
{
	CDialog::OnCancel();
}

//修改采集器状态提示信息，0:未连接,1:已连接,2:正在采集,3:正在连接，4:正在断开，其他:未知
void FingerReg::chgEgeState(int state)
{
	LPCTSTR info;
	switch (state)
	{
	case 0:info = _T("未连接"); break;
	case 1:info = _T("已连接"); break;
	case 2:info = _T("正在采集"); break;
	case 3:info = _T("正在连接"); break;
	case 4:info = _T("正在断开"); break;
	default:
		info = _T("未知");
		break;
	}
	m_textEngState.SetWindowTextW(info);
}

//获取输入的名字
//返回wchar_t长度
int FingerReg::GetInputName(wchar_t ** name)
{
	CEdit* editor;
	editor = (CEdit*)GetDlgItem(IDC_EDIT_NAME);
	//取值
	CString str;
	editor->GetWindowTextW(str);
	LPCWSTR c = str.GetString();
	int length = str.GetLength();
	*name = (wchar_t*)malloc(sizeof(wchar_t)*(length + 1));
	(*name)[length] = '\0';
	memcpy(*name, c, length * sizeof(wchar_t));
	//赋值为空
	editor->SetWindowText(_T(""));
	return length;
}



//采集器接收到图像后的响应函数的关联关系 //（在本系统中）相当于主函数
void FingerReg::OnImageReceivedZkfpengx(BOOL FAR* AImageValid)
{
	//HDC hdc = this->GetDC()->m_hDC;
	//MessageBox(_T("指纹采集器工作", _T("提示")));
	//int x = 0, y = 0;//图像绘制左上角坐标
	//long width = m_zkfpEng.get_ImageWidth();
	//long height = m_zkfpEng.get_ImageHeight();
	//m_zkfpEng.PrintImageAt(int(hdc), x, y, width, height);


	unsigned char * ucImg = NULL;  //每次处理图像的缓冲区
	unsigned char * ucDstImg = NULL; //目标图像缓冲区


	int iWidth;
	int iHeight;
	int iDepth;

	isRefresh = true;
	chgEgeState(2);
	if (*AImageValid == true)
		m_zkfpEng.SaveBitmap(FingerSolve::ToWideChar(PATHCAT(DEBUGFINGERS,"capture.bmp")));
	chgEgeState(1);
	if (imageShowRadioCheckedIndex == 1) ShowImageInCtrl1(m_picCtrl,PATHCAT(DEBUGFINGERS,"capture.bmp"));


	wchar_t * name;
	if (GetInputName(&name) == 0&& flag == 1) {//如果是登记操作，名字不能为空
		MessageBoxW(_T("警告：名字不能为空"));
		return;
	};

	//1、打开图像 初始化目标图像缓冲区=>复制图像到目标图像缓冲区
	int result = FingerSolve::readBMP2buf(PATHCAT(DEBUGFINGERS,"capture.bmp"), &ucImg, iWidth, iHeight, iDepth);

	if (result == -1)return;

	//2、中值滤波 
	ucDstImg = (byte*)malloc(sizeof(byte)*iWidth*iHeight);
	FingerSolve::MidFilter(ucImg, ucDstImg, iWidth, iHeight);
	//memcpy(ucImg, ucDstImg, sizeof(ucDstImg));
	////2次中值滤波
	//ucDstImg = (byte*)malloc(sizeof(byte)*iWidth*iHeight);
	//FingerSolve::MidFilter(ucImg, ucDstImg, iWidth, iHeight);
	if (imageShowRadioCheckedIndex == 2) {
		FingerSolve::SaveDataToImageFile(PATHCAT(DEBUGFINGERS,"capture.bmp"), PATHCAT(DEBUGFINGERS, "capture-midFileter.bmp"), ucDstImg);
		ShowImageInCtrl1(m_picCtrl, PATHCAT(DEBUGFINGERS, "capture-midFileter.bmp"));
	}
	memcpy(ucImg, ucDstImg, sizeof(ucDstImg));

	//3、均衡化
	FingerSolve::HistoNormalize(ucImg, ucDstImg, iWidth, iHeight);
	if (imageShowRadioCheckedIndex == 3)
	{
		FingerSolve::SaveDataToImageFile(PATHCAT(DEBUGFINGERS,"capture.bmp"), PATHCAT(DEBUGFINGERS, "capture-histonormalize.bmp"), ucDstImg);
		ShowImageInCtrl1(m_picCtrl, PATHCAT(DEBUGFINGERS, "capture-histonormalize.bmp"));
	}
	memcpy(ucImg, ucDstImg, sizeof(ucDstImg));

	//4、指纹图像方向计算、指纹脊线计算 3.11
	float *tmpDirections = new float[iWidth*iHeight];
	float *directions = new float[iWidth*iHeight];
	float *frequence = new float[iWidth*iHeight];
	byte *msk = new byte[iWidth*iHeight];
	byte * ucImgGaborEnhenced = new byte[iWidth*iHeight];//Gabor增强后的图像缓冲区

	FingerSolve::ImgDirection(ucImg, tmpDirections, iWidth, iHeight); //图像方向计算
	FingerSolve::DircLowPass(tmpDirections, directions, iWidth, iHeight); // 低通滤波
	if (imageShowRadioCheckedIndex == 4) {
		FingerSolve::SaveDataToImageFile(PATHCAT(DEBUGFINGERS,"capture.bmp"), PATHCAT(DEBUGFINGERS, "capture-direction.bmp"), directions, iWidth, iHeight, iDepth, 100);
		ShowImageInCtrl1(m_picCtrl, PATHCAT(DEBUGFINGERS, "capture-direction.bmp"));
	}

	//5、频率计算 3.12
	FingerSolve::Frequency(ucImg, directions, frequence, iWidth, iHeight);//指纹脊线频率计算
	if (imageShowRadioCheckedIndex == 5)
	{
		FingerSolve::SaveDataToImageFile(PATHCAT(DEBUGFINGERS,"capture.bmp"), PATHCAT(DEBUGFINGERS, "capture-frequence.bmp"), frequence, iWidth, iHeight, iDepth, 1000);
		ShowImageInCtrl1(m_picCtrl, PATHCAT(DEBUGFINGERS, "capture-frequence.bmp"));
	}

	//6 掩码计算 3.13
	FingerSolve::GetMask(ucImg, directions, frequence, msk, iWidth, iHeight);//掩码计算
	if (imageShowRadioCheckedIndex == 6)
	{
		FingerSolve::SaveDataToImageFile(PATHCAT(DEBUGFINGERS,"capture.bmp"), PATHCAT(DEBUGFINGERS, "capture-mask.bmp"), msk);
		ShowImageInCtrl1(m_picCtrl, PATHCAT(DEBUGFINGERS, "capture-mask.bmp"));
	}

	//7 Gabor增强 3.14
	FingerSolve::GaborEnhance(ucImg, directions, frequence, msk, ucImgGaborEnhenced, iWidth, iHeight);
	if (imageShowRadioCheckedIndex == 7) {
		FingerSolve::SaveDataToImageFile(PATHCAT(DEBUGFINGERS,"capture.bmp"), PATHCAT(DEBUGFINGERS, "capture-gaborEnhenced.bmp"), ucImgGaborEnhenced);
		ShowImageInCtrl1(m_picCtrl, PATHCAT(DEBUGFINGERS, "capture-gaborEnhenced.bmp"));
	}

	//8 指纹图像二值化 3.16
	byte * BinImg = ucDstImg; //直接利用现有内存区，节省内存
	byte * GrayImg = ucImg; //直接利用现有内存去，节省内存

	FingerSolve::BinaryImg(ucImgGaborEnhenced, BinImg, iWidth, iHeight, 80);
	FingerSolve::BinaryToGray(BinImg, GrayImg, iWidth, iHeight);
	if (imageShowRadioCheckedIndex == 8)
	{
		FingerSolve::SaveDataToImageFile(PATHCAT(DEBUGFINGERS,"capture.bmp"), PATHCAT(DEBUGFINGERS, "capture-bin2gray.bmp"), GrayImg);
		ShowImageInCtrl1(m_picCtrl, PATHCAT(DEBUGFINGERS, "capture-bin2gray.bmp"));
	}

	//9 细化 3.17
	byte * ThinImg = ucImg;
	byte * ThinedGray = ucDstImg;

	FingerSolve::Thinning(BinImg, ThinImg, iWidth, iHeight, 200);
	FingerSolve::BinaryToGray(ThinImg, ThinedGray, iWidth, iHeight);
	if (imageShowRadioCheckedIndex == 9) {
		FingerSolve::SaveDataToImageFile(PATHCAT(DEBUGFINGERS,"capture.bmp"), PATHCAT(DEBUGFINGERS, "capture-thined.bmp"), ThinedGray);
		ShowImageInCtrl1(m_picCtrl, PATHCAT(DEBUGFINGERS, "capture-thined.bmp"));
	}

	delete[]tmpDirections;
	delete[]directions;
	delete[]frequence;
	delete[]msk;
	delete[]ucImgGaborEnhenced;

	//10 特征提取 3.18
	byte *ucMinuImg = new byte[iWidth*iHeight];
	byte *MinuImgGray = ucDstImg;
	int MCount;
	MCount = FingerSolve::Extract(ThinImg, ucMinuImg, iWidth, iHeight);
	FingerSolve::BinaryToGray(ucMinuImg, MinuImgGray, iWidth, iHeight);
	if (imageShowRadioCheckedIndex == 10) {
		FingerSolve::SaveDataToImageFile(PATHCAT(DEBUGFINGERS,"capture.bmp"), PATHCAT(DEBUGFINGERS, "capture-features.bmp"), MinuImgGray);
		ShowImageInCtrl1(m_picCtrl, PATHCAT(DEBUGFINGERS, "capture-features.bmp"));
	}

	//11 特征过滤 3.19
	minutiae * minutiaes = new minutiae[MCount];

	FingerSolve::MinuFilter(ucMinuImg, ThinImg, minutiaes, MCount, iWidth, iHeight);
	memset(MinuImgGray, 0, iWidth*iHeight);
	for (size_t i = 0; i < MCount; i++)
	{
		MinuImgGray[(minutiaes[i].y - 1)*iWidth + (minutiaes[i].x - 1)] = 255;
	}
	if (imageShowRadioCheckedIndex == 11) {
		FingerSolve::SaveDataToImageFile(PATHCAT(DEBUGFINGERS,"capture.bmp"), PATHCAT(DEBUGFINGERS, "capture-features_f.bmp"), MinuImgGray);
		ShowImageInCtrl1(m_picCtrl, PATHCAT(DEBUGFINGERS, "capture-features_f.bmp"));
	}

	if (flag == 1)
	{
		//登记入库
		MINURECORD * minurecord = new MINURECORD[1];
		minurecord->name = name; //从输入框获取的名字
		minurecord->minu = minutiaes;
		minurecord->minuCount = MCount;
		int r = mdb->InsertNewMinu(minurecord);

		if (r < 1) {
			MessageBoxW(_T("信息录入失败"));
			return;
		}
		MessageBoxW(_T("信息录入成功"), _T("指纹录入提示"), MB_OK);
		delete[] minurecord;
	}

	/*指纹读取部分*//*
	MINURECORD * minurecord1 = new MINURECORD;
	mdb->GetMinu(1, minurecord1);
	delete[] minurecord1;*/

	/*获取全部指纹记录*/
	if (flag == 2)
	{
		MINURECORD * minurecord2 = NULL;
		int recordCount;
		mdb->GetMinu(&minurecord2, &recordCount);

		int id = -1;
		float similarity = -1;


		for (size_t i = 0; i < recordCount; i++)
		{
			int tmpsimilarity = FingerSolve::MinuSimilarity(iWidth, iHeight, minutiaes, MCount, minurecord2[i].minu, minurecord2[i].minuCount
			);
			if (tmpsimilarity > similarity) {
				similarity = tmpsimilarity;
				id = i;
			}
		}

		if (similarity < 0.4)
		{
			MessageBoxW(_T("识别失败"));
		}
		else
		{
			char buffer[200];
			//sprintf_s(buffer, "%s:%d,%s:%d", "识别精度：", similarity, "识别结果：", minurecord2[id].name);
			MessageBoxW(minurecord2[id].name, _T("识别结果"), MB_OK);
			//		MessageBoxW((wchar_t*)&similarity, _T("识别精度"), MB_OK);
		}
		delete[] minurecord2;
	}
	
	free(ucMinuImg);
	free(ucImg);
	free(ucDstImg);
	free(minutiaes);
}


//停止采集器
void FingerReg::OnBnClickedButtonStopengine()
{
	stopEngine();
}

//显示图像文件到图像控件 使用版本
int ShowImageInCtrl1(CStatic &picCtrl, char*filename)
{
	//isRefresh = false;
	CImage image;
	HRESULT hResult = image.Load(FingerSolve::ToWideChar(filename));
	if (FAILED(hResult) || image.IsNull()) {
		return -1;
	}

	picCtrl.ModifyStyle(SS_TYPEMASK, SS_BITMAP);

	CDC *pDc = picCtrl.GetWindowDC();
	SetStretchBltMode(pDc->m_hDC, STRETCH_HALFTONE);

	picCtrl.SetBitmap(image);

	//picCtrl.Invalidate(false);

	image.Destroy();
	picCtrl.ReleaseDC(pDc);

	return 0;
}

void OnPaint()
{
	if (!isRefresh)
		isShowIng = true;
}


//无采集器测试按钮处理函数
void FingerReg::OnBnClickedButtontest()
{
	flag = 1;
	BOOL FAR aImageValid = false;//保证不会调用采集器图像保存函数SaveImage()
	OnImageReceivedZkfpengx(&aImageValid);
}


//选择显示的图像的RadioButton处理函数
void FingerReg::OnBnClickedRadioShowImage()
{
	if (IsDlgButtonChecked(IDC_RADIO_SOURCE))imageShowRadioCheckedIndex = 1; //原图
	else if (IsDlgButtonChecked(IDC_RADIO2))imageShowRadioCheckedIndex = 2;//中值滤波
	else if (IsDlgButtonChecked(IDC_RADIO3))imageShowRadioCheckedIndex = 3;//均衡化
	else if (IsDlgButtonChecked(IDC_RADIO4))imageShowRadioCheckedIndex = 4;//脊线方向
	else if (IsDlgButtonChecked(IDC_RADIO5))imageShowRadioCheckedIndex = 5;//脊线频率
	else if (IsDlgButtonChecked(IDC_RADIO6))imageShowRadioCheckedIndex = 6;//掩码计算
	else if (IsDlgButtonChecked(IDC_RADIO7))imageShowRadioCheckedIndex = 7;//Gabor滤波
	else if (IsDlgButtonChecked(IDC_RADIO8))imageShowRadioCheckedIndex = 8;//二值化
	else if (IsDlgButtonChecked(IDC_RADIO9))imageShowRadioCheckedIndex = 9;//细化
	else if (IsDlgButtonChecked(IDC_RADIO10))imageShowRadioCheckedIndex = 10;//特征
	else if (IsDlgButtonChecked(IDC_RADIO11))imageShowRadioCheckedIndex = 11;//过滤特征

}


void FingerReg::OnBnClickedButtontestvalid()
{
	// TODO: 在此添加控件通知处理程序代码
	if (flag == 2)
	{
		flag = 1;
		MessageBoxW(_T("已经切换到登记模式"));
	}
	else
	{
		flag = 2;
		MessageBoxW(_T("已经切换到匹配模式"));
	}
	//BOOL FAR aImageValid = false;//保证不会调用采集器图像保存函数SaveImage()
	//OnImageReceivedZkfpengx(&aImageValid);
}


BOOL FingerReg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// TODO:  在此添加额外的初始化
	initEngine();

	return TRUE;  // return TRUE unless you set the focus to a control
				  // 异常: OCX 属性页应返回 FALSE
}
