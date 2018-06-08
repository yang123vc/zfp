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
	InitMinuDB();//��ʼ��ָ�����ݿ�
}

int db::GetNewIndexDB()
{
	return 0;
}

/*
�����µ�ָ�ƵǼ�
param:
	name ����
	minu ָ��������
	minuCount ָ����������
*/
int db::InsertNewMinu(MINURECORD *record)
{
	opendb();
	sqlite3_stmt *stat;
	sqlite3_prepare(_db, "insert into MINUTIAES(ID,NAME, MINU,MINUCOUNT) values(null,?, ? ,?)", -1, &stat, 0);
	sqlite3_bind_text16(stat, 1, record->name, sizeof(record->name), NULL);//��һ���������Ǽ���������
	sqlite3_bind_blob(stat, 2, record->minu, sizeof(minutiae)*record->minuCount, NULL);//�ڶ���������ָ��������
	sqlite3_bind_int(stat, 3, record->minuCount);//�ڶ���������ָ��������
	int r = sqlite3_step(stat);//ִ��statment
	sqlite3_finalize(stat);//�ͷ�statment
	closedb();
	return r;//����ɹ�
}

/*
�����ֻ�ȡָ�ƵǼ���Ϣ
*/
int db::GetMinu(wchar_t name, MINURECORD * record)
{
	return 0;
}

/*
��ָ����Ż�ȡָ�ƵǼ���Ϣ
return -1,����return 0 ��ȷ;
*/
int db::GetMinu(int index, MINURECORD * record)
{
	sqlite3_stmt *stmt;
	opendb();
	sqlite3_prepare(_db, "select * from MINUTIAES where ID=1;", -1, &stmt, NULL);
	//sqlite3_bind_int(stmt, 1, index);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		//ִ��stmt
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


//�����ݿ�
int db::opendb()
{
	char *zErrMsg = 0;
	int rc;

	rc = sqlite3_open(_dbname, &_db);

	if (rc) {
		fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(_db));
		return -1;//��ʧ��
	}
	else {
		fprintf(stdout, "Opened database successfully\n");
		return 0;//�򿪳ɹ�
	}

}

//�ر����ݿ�
int db::closedb()
{
	sqlite3_close(_db);
	return 0;
}

/*
��ʼ��ָ�����ݿ⣻
����ָ�����ݿ�
*/

int db::InitMinuDB()
{

	char *zErrMsg = 0;
	int  rc;
	char *sql;

	if (opendb() != 0) { return -1; }//�����ݿ�ʧ��

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