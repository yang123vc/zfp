
// HFingerPDlg.cpp: 实现文件
//

#include "stdafx.h"
#include "HFingerP.h"
#include "HFingerPDlg.h"
#include "afxdialogex.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CHFingerPDlg 对话框



CHFingerPDlg::CHFingerPDlg(CWnd* pParent /*=nullptr*/)
	: CDialogEx(IDD_HFINGERP_DIALOG, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CHFingerPDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);

}

BEGIN_MESSAGE_MAP(CHFingerPDlg, CDialogEx)
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDCANCEL, &CHFingerPDlg::OnBnClickedCancel)
//	ON_BN_CLICKED(IDC_BUTTON2, &CHFingerPDlg::OnBnClickedButton2)
ON_BN_CLICKED(IDC_BUTTON2, &CHFingerPDlg::OnBnClickedButton2)
ON_BN_CLICKED(IDC_BUTTON1, &CHFingerPDlg::OnBnClickedButton1)
END_MESSAGE_MAP()


// CHFingerPDlg 消息处理程序

BOOL CHFingerPDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// 设置此对话框的图标。  当应用程序主窗口不是对话框时，框架将自动
	//  执行此操作
	SetIcon(m_hIcon, TRUE);			// 设置大图标
	SetIcon(m_hIcon, FALSE);		// 设置小图标

	// TODO: 在此添加额外的初始化代码
	


	return TRUE;  // 除非将焦点设置到控件，否则返回 TRUE
}

// 如果向对话框添加最小化按钮，则需要下面的代码
//  来绘制该图标。  对于使用文档/视图模型的 MFC 应用程序，
//  这将由框架自动完成。

void CHFingerPDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // 用于绘制的设备上下文

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// 使图标在工作区矩形中居中
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// 绘制图标
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

//当用户拖动最小化窗口时系统调用此函数取得光标
//显示。
HCURSOR CHFingerPDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}



void CHFingerPDlg::OnBnClickedCancel()
{
	//退出按钮处理事件
	CDialogEx::OnCancel();
}



void CHFingerPDlg::OnBnClickedButton2() //指纹登记按钮处理事件
{
	// 指纹登记函数，点击按钮“指纹登记”弹出指纹登记对话框，处理指纹登记。
	FingerReg dlg_fingerreg; //登记对话框对象
	INT_PTR nResponse = dlg_fingerreg.DoModal();
	if (nResponse == IDOK)
	{
		// TODO: 在此放置处理何时用
		//  “确定”来关闭对话框的代码
	}
	else if (nResponse == IDCANCEL)
	{
		// TODO: 在此放置处理何时用
		//  “取消”来关闭对话框的代码
	}
	else if (nResponse == -1)
	{
		TRACE(traceAppMsg, 0, "警告: 对话框创建失败，应用程序将意外终止。\n");
		TRACE(traceAppMsg, 0, "警告: 如果您在对话框上使用 MFC 控件，则无法 #define _AFX_NO_MFC_CONTROLS_IN_DIALOGS。\n");
	}
	
}


void CHFingerPDlg::OnBnClickedButton1()
{
	// TODO: 在此添加控件通知处理程序代码
	// 指纹登记函数，点击按钮“指纹登记”弹出指纹登记对话框，处理指纹登记。
	//FingerReg1 dlg_fingerreg; //登记对话框对象
	//INT_PTR nResponse = dlg_fingerreg.DoModal();
	//if (nResponse == IDOK)
	//{
	//	// TODO: 在此放置处理何时用
	//	//  “确定”来关闭对话框的代码
	//}
	//else if (nResponse == IDCANCEL)
	//{
	//	// TODO: 在此放置处理何时用
	//	//  “取消”来关闭对话框的代码
	//}
	//else if (nResponse == -1)
	//{
	//	TRACE(traceAppMsg, 0, "警告: 对话框创建失败，应用程序将意外终止。\n");
	//	TRACE(traceAppMsg, 0, "警告: 如果您在对话框上使用 MFC 控件，则无法 #define _AFX_NO_MFC_CONTROLS_IN_DIALOGS。\n");
	//}
}
