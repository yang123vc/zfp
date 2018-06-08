#pragma once
#include "sqlite3.h"
#include "FingerSolve.h"

struct MINURECORD {
	int ID;
	wchar_t* name;
	struct minutiae *minu;
	int minuCount;
};

class db
{
public:
	db();
	db(const char * dbname);
	int GetNewIndexDB();
	int InsertNewMinu(MINURECORD *record);
	int GetMinu(wchar_t name, MINURECORD * record);
	int GetMinu(MINURECORD ** record, int * recordCount);
	int GetMinu(int index, MINURECORD * record);
	~db();
private:
	sqlite3 * _db;
	char * _dbname;
	int opendb();
	int closedb();
	int InitMinuDB();
};

