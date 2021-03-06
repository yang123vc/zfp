#pragma once
#include "CZKFPEngX.h"
#include "FingerSolve.h"
#include "db.h"
#include <string>

// FingerReg 对话框

class FingerReg : public CDialog
{
	DECLARE_DYNAMIC(FingerReg)

public:
	CZKFPEngX m_zkfpEng;//指纹采集器对象
	CStatic m_picCtrl;
	CButton m_btnInitEngine,m_btnStopEngine;
	CStatic m_textEngState;

	FingerReg(CWnd* pParent = nullptr);   // 标准构造函数
	virtual ~FingerReg();


private:
	int imageShowRadioCheckedIndex;
	db *mdb;

// 对话框数据
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_FINGERREG_DIALOG };
#endif

protected:
	int flag; //1登记，2识别
	void initEngine();
	void stopEngine();
	void chgEgeState(int);
	int GetInputName(wchar_t ** name);//获取输入的名字
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持
	DECLARE_EVENTSINK_MAP()
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedCancel();
	afx_msg void OnBnClickedButtonInitengenine();
	afx_msg void OnBnClickedButtonStopengine();

	afx_msg void OnImageReceivedZkfpengx(BOOL FAR* AImageValid);
	afx_msg void OnBnClickedButtontest();
	afx_msg void OnBnClickedRadioShowImage();

	afx_msg void OnBnClickedButtontestvalid();
	virtual BOOL OnInitDialog();
};
