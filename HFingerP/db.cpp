#include "stdafx.h"
#include "db.h"


db::db()
{

}

db::db(const char * dbname)
{
	_dbname = (char*)dbname;
	opendb();
	closedb();
	InitMinuDB();//初始化指纹数据库
}

int db::GetNewIndexDB()
{
	return 0;
}

/*
插入新的指纹登记
param:
	name 名字
	minu 指纹特征点
	minuCount 指纹特征点数
*/
int db::InsertNewMinu(MINURECORD *record)
{
	opendb();
	sqlite3_stmt *stat;
	sqlite3_prepare(_db, "insert into MINUTIAES(ID,NAME, MINU,MINUCOUNT) values(null,?, ? ,?)", -1, &stat, 0);
	sqlite3_bind_text16(stat, 1, record->name, sizeof(record->name), NULL);//第一个参数（登记人姓名）
	sqlite3_bind_blob(stat, 2, record->minu, sizeof(minutiae)*record->minuCount, NULL);//第二个参数（指纹特征）
	sqlite3_bind_int(stat, 3, record->minuCount);//第二个参数（指纹特征）
	int r = sqlite3_step(stat);//执行statment
	sqlite3_finalize(stat);//释放statment
	closedb();
	return r;//插入成功
}

/*
按名字获取指纹登记信息
*/
int db::GetMinu(wchar_t name, MINURECORD * record)
{
	return 0;
}

/*
按指纹序号获取指纹登记信息
return -1,错误，return 0 正确;
*/
int db::GetMinu(int index, MINURECORD * record)
{
	sqlite3_stmt *stmt;
	opendb();
	sqlite3_prepare(_db, "select * from MINUTIAES where ID=1;", -1, &stmt, NULL);
	//sqlite3_bind_int(stmt, 1, index);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		//执行stmt
		record->ID = sqlite3_column_int(stmt, 0);
		int nameSize = sizeof((wchar_t*)sqlite3_column_text16(stmt, 1))/2;
		record->name = new wchar_t [nameSize+2];
		record->name[nameSize] = '\0';
		memcpy(record->name,(wchar_t*)sqlite3_column_text16(stmt,1),nameSize*2);
		record->minuCount = sqlite3_column_int(stmt, 3);
		record->minu = new minutiae[record->minuCount];
		memcpy(record->minu, (minutiae*)sqlite3_column_blob(stmt, 2), record->minuCount*sizeof(minutiae));
		break;
	}
	sqlite3_finalize(stmt);
	closedb();
	return 0 ;
}


//打开数据库
int db::opendb()
{
	char *zErrMsg = 0;
	int rc;

	rc = sqlite3_open(_dbname, &_db);

	if (rc) {
		fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(_db));
		return -1;//打开失败
	}
	else {
		fprintf(stdout, "Opened database successfully\n");
		return 0;//打开成功
	}

}

//关闭数据库
int db::closedb()
{
	sqlite3_close(_db);
	return 0;
}

/*
初始化指纹数据库；
创建指纹数据库
*/

int db::InitMinuDB()
{

	char *zErrMsg = 0;
	int  rc;
	char *sql;

	if (opendb() != 0) { return -1; }//打开数据库失败

	/* Create SQL statement */
	sql = "CREATE TABLE MINUTIAES(ID INTEGER PRIMARY KEY AUTOINCREMENT,NAME TEXT NOT NULL,MINU BLOB NOT NULL,MINUCOUNT INT NOT NULL);";

	/* Execute SQL statement */
	rc = sqlite3_exec(_db, sql, NULL, 0, &zErrMsg);
	if (rc != SQLITE_OK) {
		fprintf(stderr, "SQL error: %s\n", zErrMsg);
		sqlite3_free(zErrMsg);
	}
	else {
		fprintf(stdout, "Table created successfully\n");
	}

	closedb();

	return 0;
}

db::~db()
{
	sqlite3_close(_db);
}