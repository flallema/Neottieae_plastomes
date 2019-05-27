//============================================================================
// Name        : SearchRedondance.cpp
// Author      : Maria Tchoumakov
// Version     : 1.7.6
// Description : Identification of redondant sequences from rawdata by couple Read1-Read2or in one single read.
//
// UPDATE (by Sivasangari Nandy) - 07/03/2016 : new version. Insertion in DB (var\tval)
//
// INPUT : DBData.tmp, Connector.tmp,ID_Processus 
// OUTPUT : uniq sequence files, redundancy file with the frequency of each redundant sequence
//============================================================================

using namespace std;

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <istream>
#include <string.h>
#include <fstream>
#include <ostream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <map>
#include <list>
#include <exception>
#include <mysql.h>
#include <time.h>
map<string,string> ReadDataFile(string);
string FormatDate();
void createRepeatHashtable(map<string,string>,string,string,map<string,string>,string);
void createRepeatHashtablepairend(map<string,string>,string,string,map<string,string>,string);

class Test{
    public:
       int calcPhredScore(string);
       int calcPhredScoreSanger(string);//CH
       int calcPhredScoreIndex(string,int);
       int calcPhredScoreIndexSanger(string,int); //CH
       float calcPhredScoreMoyen(string,int, int);
       float calcPhredScoreMoyenSanger(string,int, int); //CH
       string extractString(string,int, int);
       int num;
       string header1;
       string seq;
       string header2;
       string quality;
       int score;
       int indice;
       
       Test (string laseq)
       {
           seq=laseq;
       }
       
       Test (int Num, string Header1, string Seq, string Header2, string Quality, int Indice){
           num = Num;
           header1 = Header1;
           seq=Seq;
           header2 = Header2;
           quality = Quality;
           score=this->
           calcPhredScore(Quality);
           indice=Indice;
       }
};

int Test::calcPhredScore(string scoreline){
    int score1=0;
    stringstream ss;
    char ascii;
    int val;
    ss<<scoreline;
    while(ss.good()){
        ss>>ascii;
        val=(int)ascii;
        score1=score1+val-64;
    }
    return score1;
}

int Test::calcPhredScoreSanger(string scoreline ) 
{
    int score1=0;
    stringstream ss;
    char ascii;
    int val;
    ss<<scoreline;
    while(ss.good()){
        ss>>ascii;
        val=(int)ascii;
        score1=score1+val-33;
    }
    return score1;
}

int Test::calcPhredScoreIndex(string seq, int index){

	int size = seq.size() + 1;
	char * scoretab = new char[ size ];
	strncpy( scoretab, seq.c_str(), size);
	int score = (int) scoretab[index]-64;
	delete [] scoretab;
	return score;
}

int Test::calcPhredScoreIndexSanger(string seq, int index){

	int size = seq.size() + 1;
	char * scoretab = new char[ size ];
	strncpy( scoretab, seq.c_str(), size);
	int score = (int) scoretab[index]-33;
	delete [] scoretab;
	return score;
}
float Test::calcPhredScoreMoyen(string seq, int debut, int fin){

	int size = seq.size() + 1;
	char * scoretab = new char[ size ];
	strncpy( scoretab, seq.c_str(), size);
	int somme=0;
	float moyenne=0;
	for (int i=debut; i<=fin; ++i){
		int score = (int) scoretab[i]-64;
		somme = somme + score;
	}
	if (fin-debut+1>0){
		moyenne = somme / (fin-debut+1);
	}
	else {
		cout << "erreur moyenne, div0" << endl;
		exit (0);
	}

	delete [] scoretab;

	return moyenne;
}

float Test::calcPhredScoreMoyenSanger(string seq, int debut, int fin){

	int size = seq.size() + 1;
	char * scoretab = new char[ size ];
	strncpy( scoretab, seq.c_str(), size);
	int somme=0;
	float moyenne=0;
	for (int i=debut; i<=fin; ++i){
		int score = (int) scoretab[i]-33;
		somme = somme + score;
	}
	if (fin-debut+1>0){
		moyenne = somme / (fin-debut+1);
	}
	else {
		cout << "erreur moyenne, div0" << endl;
		exit (0);
	}

	delete [] scoretab;

	return moyenne;
}

string Test::extractString(string seq, int debut, int fin){
	int size = seq.size() + 1;
		char * scoretab = new char[ size ];
		strncpy( scoretab, seq.c_str(), size);
		string newstring="";
		for (int i=debut; i<=fin; ++i){
			newstring+=scoretab[i];
		}
		delete [] scoretab;

		return newstring;
}

bool compare_readrep (string first, string second)
{
  stringstream  ss1;
  stringstream  ss2;
  ss1<<first;
  ss2<<second;

  int a, b;
  ss1>>a;
  ss2>>b;
  if ((ss1.good())&&(ss2.good())){
        if (a>b) return true;
        else return false;
  }
  else return false;
}

int main(int argc,char* argv[]){
	string Id_Processus=argv[1];
	string commonName=argv[2];
	string FileName=argv[3];
	
	string FileConnectTmp=argv[4];  
	map<string,string> SampleData=ReadDataFile(FileName);
	map<string,string>Connector=ReadDataFile(FileConnectTmp);
	if(((SampleData["PairEnd"]).compare("PE"))==0){
		createRepeatHashtablepairend(SampleData,commonName,FileName,Connector,Id_Processus);
	}
	else if(((SampleData["PairEnd"]).compare("PE"))!=0){
		createRepeatHashtable( SampleData,commonName,FileName,Connector,Id_Processus);	
	}
	return 0;
}


string FormatDate(){
	time_t rawtime;
	time(&rawtime);
	struct tm * ptm;
	//2010-03-06 02:00:00
	//timestamp=localtime(&rawtime);

	ptm = gmtime ( &rawtime );
	int mon=(ptm->tm_mon)+1;
	stringstream ss;
	ss<<mon;
	string month=ss.str();
	if(month.size()==1){
		month="0"+month;
	}
	stringstream si;
	si<<((ptm->tm_year)+1900);
	string year=si.str();
	stringstream sm;
	sm<<ptm->tm_mday;
	string day=sm.str();
	if(day.size()==1){
		day="0"+day;
	}
	stringstream sh;
	sh<<((ptm->tm_hour)+2);
	string hour=sh.str();
	if(hour.size()==1){
		hour="0"+hour;
	}
	stringstream smin;
	smin<<ptm->tm_min;
	string min=smin.str();
	if(min.size()==1){
		min="0"+min;
	}
	stringstream ssec;
	ssec<<ptm->tm_sec;
	string sec=ssec.str();
	if(sec.size()==1){
		sec="0"+sec;
	}
	string timestr=year+":"+month+":"+day+" "+hour+":"+min+":"+sec;
	return timestr;
}
void Split(vector<string>& vecteur, string chaine, const char * separateur)
{
	vecteur.clear();
	string::size_type stTemp = chaine.find(separateur);
	while(stTemp != string::npos)
	{
		vecteur.push_back(chaine.substr(0, stTemp));
		chaine = chaine.substr(stTemp + 1);
		stTemp = chaine.find(separateur);
	}
}
map<string,string>ReadDataFile(string FileName){
	string line;
	ifstream fin;
	vector<string> words;
	map<string,string>SampleData;
	map<string,string>::iterator i;
	pair<map<string,string>::iterator,bool> ret;
	fin.open(FileName.c_str());
	if(fin.is_open()){

	}
	else{
		cout<<1<<endl;
		cout<<"probleme à l'ouverture du fichier "<<FileName<<endl;
		exit(1);
	}
	int count=0;
	while(!fin.eof()){
		getline(fin,line);
		count++;
		Split(words,line," ");
		SampleData.insert( pair<string,string>(words[0],words[1]) );
	}
	fin.close();

	return SampleData;
}
void createRepeatHashtablepairend(map<string,string>SampleData,string commonName,string FileDB,map<string,string>Connector,string ID_Processus){
    // var for DB
	string step="";
	string timestamp="";
	string variable="";
	string value="";
	string type="";
	string comment="";
	string line="";
	string step_number="";
	string query="";
	int linenum=0;	
	string ErrorCode="";
	string ErrorMessage="";
	MYSQL *connect;
	MYSQL *status;
	MYSQL_ROW row;
	MYSQL_RES *result;
	string readname="";
	string server=Connector["host"];
	string database=Connector["dbname"];
	string user=Connector["user"];
	string pass=Connector["password"];
	string timestamp_deb=FormatDate();
	connect=mysql_init(NULL);
	char* cmd;
	ifstream fin1;
    ifstream fin2;
    ofstream fout1;
    ofstream fout2;
    ofstream fout3;
	ofstream fout5;
    int comptread=0;
	list<int>listindex;
    map<string,int> maplines;
    map<string,int> scorelines;
    map<string,int> indicelines;
    list<string> listlines;
	Test* membre = new Test;
	string Read1=SampleData["R1_Filename"];
	string Read2=SampleData["R2_Filename"];	
	string AnalyseREP=SampleData["REPANAL"];
	string DecomREP=SampleData["REPDECOMP"];	
	string CasavaVersion=SampleData["Casava_Version"];
	string Fileout1="uniq_"+Read1;
	string Fileout2="uniq_"+Read2;
	string Fileout3="redondances_"+commonName+".txt";
	string Path1=DecomREP+"/"+Read1;
	string Path2=DecomREP+"/"+Read2;

	// INSERTION INTO DB
	string queryStep="call PC_GetNextStep(@idProcessus:='"+ID_Processus+"');";

	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,queryStep.c_str())){
			fprintf(stderr, "%s\n", mysql_error(connect));
			cout<<cmd<<endl;
			ErrorMessage="probleme de connexion à la Database PC_GetNextStep";
			ErrorCode="1";
			exit(1);
		}
		result=mysql_store_result(connect);
              
		if((row=mysql_fetch_row(result))){
			step_number=row[0];
		}
		else{
			cout<<"pas de resultats de requete"<<endl;
		}
	}
	else{
		ErrorMessage="probleme de connexion à la Database PC_GetNextStep";
		ErrorCode="1";
		cout<<ErrorCode<<endl;
		cout<<ErrorMessage<<endl;
		exit(1);
	}	
	cout<<"step num\t"<< step_number<<endl;
	cout<<"DEBUT DU STEP\t"<< timestamp_deb<<endl;
	
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream sif;
	sif<<linenum;	
	mysql_free_result(result);
	step="SearchRedondance";
	timestamp=FormatDate();
	variable="TIME_START";
	value=timestamp_deb;
	type="DATETYPE";
	comment="Temps du debut du step";
	line=sif.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))	{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			ErrorMessage="probleme de connexion à la Database PC_InsertIntoDBLog";
			ErrorCode="1";			
			exit(1);
		}
		else{
			ErrorCode="0";
		}
	}
	//END INSERT INTO DB
	
	string value_checkFOUT2="";
	string value_checkFOUT1="";
	string value_checkFOUT3="";
	string value_checkR1="";
	string value_checkR2="";
	
	fout1.open((AnalyseREP+"/"+Fileout1).c_str());
	if(fout1.is_open()){
	    cout<<"Fichier "<<Fileout1<<" OK\n";
		value_checkFOUT1="Fichier "+Fileout1+" OK";
	}
	else{
		cout<<"Fichier "<<Fileout1<<" non existant\n";
		value_checkFOUT1="Fichier "+Fileout1+" non existant";
	}
	fout2.open((AnalyseREP+"/"+Fileout2).c_str());
	if(fout2.is_open()){
	    cout<<"Fichier "<<Fileout2<<" OK\n";
		value_checkFOUT2="Fichier "+Fileout2+" OK";
	}
	else{
		cout<<"Fichier"<<Fileout2<<"non existant\n";
		value_checkFOUT2="Fichier "+Fileout2+" non existant";
	}
	fout3.open((AnalyseREP+"/"+Fileout3).c_str());
	if(fout3.is_open()){
	    cout<<"Fichier "<<Fileout3<<" OK\n";
		value_checkFOUT3="Fichier "+Fileout3+" OK";
	}
	else{
		cout<<"Fichier "<<Fileout3<<" non existant\n";
		value_checkFOUT3="Fichier "+Fileout3+" non existant";
	}
	
	// check fastq sequence files
	fin1.open(Path1.c_str());
	if(fin1.is_open()){
		cout<<" Fichier "<<Read1<<" ouvert\n";
		value_checkR1=" Fichier "+Read1+" ouvert";
	}
	else{
		cout<<"probleme à l'ouverture du fichier "<<Read1<<"\n";
		value_checkR1="probleme à l'ouverture du fichier "+Read1;
	}

	fin2.open(Path2.c_str());
	if(fin2.is_open()){
		cout<<" Fichier "<<Read2<<" ouvert\n";
		value_checkR2=" Fichier "+Read2+" ouvert";
	}
	else{
		cout<<"probleme à l'ouverture du fichier "<<Read2<<"\n";
		value_checkR2="probleme à l'ouverture du fichier "+Read2;
	}
	
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream cr1;
	cr1<<linenum;	
	timestamp=FormatDate();
	variable="Check_FILE";
	value=value_checkR1;
	type="TEXT";
	comment="controle ouverture de fichier";
	line=cr1.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))	{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			ErrorMessage="probleme de connexion à la Database PC_InsertIntoDBLog";
			ErrorCode="1";			
			exit(1);
		}
		else{
			ErrorCode="0";
		}
	}
	//END INSERT INTO DB	
	
	linenum++;
	stringstream cr2;
	cr2<<linenum;	
	timestamp=FormatDate();
	variable="Check_FILE";
	value=value_checkR2;
	type="TEXT";
	comment="controle ouverture de fichier";
	line=cr2.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))	{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			ErrorMessage="probleme de connexion à la Database PC_InsertIntoDBLog";
			ErrorCode="1";			
			exit(1);
		}
		else{
			ErrorCode="0";
		}
	}
	//END INSERT INTO DB	

	linenum++;
	stringstream checkF1;
	checkF1<<linenum;	
	timestamp=FormatDate();
	variable="Open_FILE";
	value=value_checkFOUT1;
	type="TEXT";
	comment="ouverture fichier "+Fileout1;
	line=checkF1.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))	{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			ErrorMessage="probleme de connexion à la Database PC_InsertIntoDBLog";
			ErrorCode="1";			
			exit(1);
		}
		else{
			ErrorCode="0";
		}
	}
	//END INSERT INTO DB	
	
	linenum++;
	stringstream checkF2;
	checkF2<<linenum;	
	timestamp=FormatDate();
	variable="Open_FILE";
	value=value_checkFOUT2;
	type="TEXT";
	comment="ouverture fichier "+Fileout2;
	line=checkF2.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))	{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			ErrorMessage="probleme de connexion à la Database PC_InsertIntoDBLog";
			ErrorCode="1";			
			exit(1);
		}
		else{
			ErrorCode="0";
		}
	}
	//END INSERT INTO DB
	
	linenum++;
	stringstream checkF3;
	checkF3<<linenum;	
	timestamp=FormatDate();
	variable="Open_FILE";
	value=value_checkFOUT3;
	type="TEXT";
	comment="ouverture fichier "+Fileout3;
	line=checkF3.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))	{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			ErrorMessage="probleme de connexion à la Database PC_InsertIntoDBLog";
			ErrorCode="1";			
			exit(1);
		}
		else{
			ErrorCode="0";
		}
	}
	//END INSERT INTO DB	
	
	//Reading of couple Read1-Read2 files
	while(fin1.good()&&fin2.good()){
		comptread++;
		string r1line1="";
		string r1line2="";
		string r1line3="";
		string r1line4="";
		string r2line1="";
		string r2line2="";
		string r2line3="";
		string r2line4="";

		getline(fin1,r1line1);
		getline(fin1,r1line2);
		getline(fin1,r1line3);
		getline(fin1,r1line4);
		getline(fin2,r2line1);
		getline(fin2,r2line2);
		getline(fin2,r2line3);
		getline(fin2,r2line4);

		if ((!r1line1.empty())&&(r1line1.substr(0,r1line1.size()-1) == r2line1.substr(0,r2line1.size()-1)))
		{
			stringstream ss1;
			stringstream ss2;
			stringstream ss3;
			stringstream ss4;
			ss1<<r1line1<<" "<<r2line1;
			ss2<<r1line2<<" "<<r2line2;
			ss3<<r1line3<<" "<<r2line3;
			ss4<<r1line4<<" "<<r2line4;

			string header1 = ss1.str();
			string seq = ss2.str();
			string header2 = ss3.str();
			string quality = ss4.str();

			if (maplines.count(seq)>0){
				maplines[seq]++;
				int score = 0;
				if(CasavaVersion.compare("1.7")==0)
				{
					score = membre->calcPhredScore(quality);
				}
				if(CasavaVersion.compare("1.8")==0)
				{
					score = membre->calcPhredScoreSanger(quality);
				}
				//selection de redondances avec meilleur score Phred,ajout de l'indice de la meilleure seq dans le hash
				if (scorelines.count(seq)>0){
					if (score > scorelines[seq]){
						scorelines[seq] = score;
						indicelines[seq] = comptread;
					}
				}
				else {
					scorelines[seq] = score;
					indicelines[seq] = comptread;
				}
			}
			else{
				maplines[seq]=1;
			}
		}
	}
	//Creation of redundancy list with the amount of the redundant sequence 
	for (map<string,int>::iterator p=maplines.begin(); p!=maplines.end();++p){
		if ((*p).second>=2){
			stringstream ss;
			ss<<(*p).second;
			ss<<" ";
			ss<<(*p).first;
			string liststring = ss.str();
			listlines.push_back(liststring);
		}
	}
	listlines.sort(compare_readrep);
	//Writing redundancy file
	cout << "Ecriture table de hash";
	for (list<string>::iterator p=listlines.begin(); p!=listlines.end();++p){
		fout3 << (*p) << endl;
	}
	cout << " : Ok"<<endl;

	fout3.close();
	maplines.clear();
	listlines.clear();
	fin1.clear();
	fin1.seekg(0);
	fin2.clear();
	fin2.seekg(0);
	
	//Wrinting of uniq read files with the best sequence score
	comptread=0;
	int countlines=0;
	while(fin1.good()&&fin2.good())
	{
		comptread++;
		stringstream ss1;
		stringstream ss2;
		stringstream ss3;
		stringstream ss4;
		string r1line1="";
		string r1line2="";
		string r1line3="";
		string r1line4="";
		string r2line1="";
		string r2line2="";
		string r2line3="";
		string r2line4="";

		getline(fin1,r1line1);
		getline(fin1,r1line2);
		getline(fin1,r1line3);
		getline(fin1,r1line4);
		getline(fin2,r2line1);
		getline(fin2,r2line2);
		getline(fin2,r2line3);
		getline(fin2,r2line4);
		
		if(!fin1.eof() && !fin2.eof()){
			ss2<<r1line2<<" "<<r2line2;
			string seq = ss2.str();
			if(indicelines.count(seq)>0 &&(comptread==indicelines[seq])){
				fout1<<r1line1<<endl;
				fout1<<r1line2<<endl;
				fout1<<r1line3<<endl;
				fout1<<r1line4<<endl;

				fout2<<r2line1<<endl;
				fout2<<r2line2<<endl;
				fout2<<r2line3<<endl;
				fout2<<r2line4<<endl;
				countlines++;
			}
			if(indicelines.count(seq)==0 ){
				fout1<<r1line1<<endl;
				fout1<<r1line2<<endl;
				fout1<<r1line3<<endl;
				fout1<<r1line4<<endl;
				
				fout2<<r2line1<<endl;
				fout2<<r2line2<<endl;
				fout2<<r2line3<<endl;
				fout2<<r2line4<<endl;
				countlines++;
			}
		}
	}
	SampleData["SEQUENCE1_ANAL"]=Fileout1;
	SampleData["SEQUENCE2_ANAL"]=Fileout2;

	cout << "PE nb sequences uniques"<<"\t"<<countlines<<endl;

	linenum++;
	stringstream seq1u;
	seq1u<<linenum;	
	timestamp=FormatDate();
	variable="SEQUENCE1_ANAL";
	value=Fileout1;
	type="TEXT";
	comment="nom de sequence nettoyee";
	line=seq1u.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))	{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			ErrorMessage="probleme de connexion à la Database PC_InsertIntoDBLog";
			ErrorCode="1";			
			exit(1);
		}
		else{
			ErrorCode="0";
		}
	}
	//END INSERT INTO DB
	
	linenum++;
	stringstream seq2u;
	seq2u<<linenum;	
	timestamp=FormatDate();
	variable="SEQUENCE2_ANAL";
	value=Fileout2;
	type="TEXT";
	comment="nom de sequence nettoyee";
	line=seq2u.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))	{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			ErrorMessage="probleme de connexion à la Database PC_InsertIntoDBLog";
			ErrorCode="1";			
			exit(1);
		}
		else{
			ErrorCode="0";
		}
	}
	//END INSERT INTO DB

	char nbUnique[20];
	sprintf (nbUnique, "%d", countlines);
	
	linenum++;
	stringstream nbu;
	nbu<<linenum;	
	timestamp=FormatDate();
	variable="NBunique";
	value=nbUnique;
	type="INTEGER";
	comment="nombre de sequences uniques";
	line=nbu.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))	{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			ErrorMessage="probleme de connexion à la Database PC_InsertIntoDBLog";
			ErrorCode="1";			
			exit(1);
		}
		else{
			ErrorCode="0";
		}
	}
	//END INSERT INTO DB

	fin1.clear();
	fin2.clear();
	fin1.close();
	fin2.close();
	fout2.close();
	fout1.close();
	
	//Wrinting of DBData file
	fout5.open(FileDB.c_str());
	if(fout5.is_open()){
		ErrorMessage="Success";
		ErrorCode="0";
	}
	else{
		ErrorMessage="probleme à l'ouverture du fichier DBData.tmp " ;
		ErrorCode=1;
		cout<<ErrorCode<<endl;
		cout<<ErrorMessage<<endl;
		exit(1);
	}
	for(map<string,string>::iterator i=SampleData.begin();i!= SampleData.end();i++){
		fout5<<(*i).first<<" "<<(*i).second<<" "<<endl;
	}
	fout5.close();	
	
	
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream endT;
	endT<<linenum;	
	string timestamp_end=FormatDate();
	cout<<"FIN DU STEP\t"<< timestamp_end<<endl;

	variable="TIME_END";
	value=timestamp_end;
	type="DATETYPE";
	comment="Temps de fin de step";
	line=endT.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";

	
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))
		{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			ErrorMessage="probleme de connexion à la Database PC_InsertIntoDBLog";
			ErrorCode="1";
			exit(1);
		}
		else{

			ErrorMessage="Success connexion DB PC_InsertIntoDBLog";
			ErrorCode="0";
		}
	}	
	mysql_close(connect);
}
void createRepeatHashtable (map<string,string>SampleData,string commonName,string FileDB,map<string,string>Connector,string ID_Processus){
    // var for DB
	string step="";
	string timestamp="";
	string variable="";
	string value="";
	string type="";
	string comment="";
	string line="";
	string step_number="";
	string query="";
	int linenum=0;	
	string ErrorCode="";
	string ErrorMessage="";
	MYSQL *connect;
	MYSQL *status;
	MYSQL_ROW row;
	MYSQL_RES *result;
	string readname="";
	string server=Connector["host"];
	string database=Connector["dbname"];
	string user=Connector["user"];
	string pass=Connector["password"];
	string timestamp_deb=FormatDate();
	cout<<"DEBUT DU STEP\t"<< timestamp_deb<<endl;
	connect=mysql_init(NULL);
	char* cmd; 
		
	ifstream fin;
    ofstream fout;
    ofstream fout2;
	ofstream fout5;	
    int comptread=0;
    map<string,int> maplines;
    map<string,int>scorelines;
    map <string,int>indicelines;
    list<string> listlines;
    Test* membre =new Test;
    list<int>listindex;
	string Read=SampleData["R1_Filename"]; 
	string DecomREP=SampleData["REPDECOMP"];	
	string Fileout1="uniq_"+Read;
	string CasavaVersion=SampleData["Casava_Version"];
	string AnalyseREP=SampleData["REPANAL"];	
	string Fileout2="redondances_"+commonName+".txt";

	// INSERTION DANS BD
	string queryStep="call PC_GetNextStep(@idProcessus:='"+ID_Processus+"');";

	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,queryStep.c_str()))
		{
			exit(1);
		}
		result=mysql_store_result(connect);
              
		if((row=mysql_fetch_row(result))){
			cout<<row[0]<<endl;
			step_number=row[0];
		}
		else{
			cout<<"pas de resultats de requete"<<endl;
		}
	}
	else{
		ErrorMessage="probleme de connexion à la Database PC_GetNextStep";
		ErrorCode="1";
		cout<<ErrorCode<<endl;
		cout<<ErrorMessage<<endl;
		exit(1);
	}		
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream sif;
	sif<<linenum;	
	mysql_free_result(result);
	step="SearchRedondance";
	timestamp=FormatDate();
	variable="TIME_START";
	value=timestamp_deb;
	type="DATETYPE";
	comment="Temps du debut du step";
	line=sif.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))
		{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			ErrorMessage="probleme de connexion à la Database PC_InsertIntoDBLog";
			ErrorCode="1";
			exit(1);
		}
		else{
			ErrorMessage="Success connexion DB PC_InsertIntoDBLog";
			ErrorCode="0";
		}
	}
	//END INSERT INTO DB	
	string value_checkFOUT2="";
	string value_checkFOUT1="";
	string value_checkR1="";

	fout.open((AnalyseREP+"/"+Fileout1).c_str());
	if(fout.is_open()){
	    cout<<"Fichier "<<Fileout1<<" OK\n";
		value_checkFOUT1="Fichier "+Fileout1+" OK";	
	}
	else{
		cout<<"Fichier "<<Fileout1<<" non existant\n";
		value_checkFOUT1="Fichier "+Fileout1+" OK";	
	}
	fout2.open((AnalyseREP+"/"+Fileout2).c_str());
	if(fout2.is_open()){
	    cout<<"Fichier "<<Fileout2<<" OK\n";
		value_checkFOUT2="Fichier "+Fileout2+" OK";		
	}
	else{
		cout<<"Fichier "<<Fileout2<<" non existant\n";
		value_checkFOUT2="Fichier "+Fileout2+" non existant";	
	}

	fin.open((DecomREP+"/"+Read).c_str());
	if(fin.is_open()){
		cout<<" Fichier "<<Read<<" ouvert\n";
		value_checkR1=" Fichier "+Read+" ouvert";		
	}
	else{
		cout<<"probleme à l'ouverture du fichier "<<Read<<"\n";
		value_checkR1="probleme à l'ouverture du fichier "+Read;		
	}
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream cr1;
	cr1<<linenum;	
	timestamp=FormatDate();
	variable="Check_FILE";
	value=value_checkR1;
	type="TEXT";
	comment="controle ouverture de fichier";
	line=cr1.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))
		{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			ErrorMessage="probleme de connexion à la Database PC_InsertIntoDBLog";
			ErrorCode="1";
			exit(1);
		}
		else{
			ErrorMessage="Success connexion DB PC_InsertIntoDBLog";
			ErrorCode="0";
		}
	}
	//END INSERT INTO DB	
	linenum++;
	stringstream checkF1;
	checkF1<<linenum;	
	timestamp=FormatDate();
	variable="Open_FILE";
	value=value_checkFOUT1;
	type="TEXT";
	comment="ouverture fichier "+Fileout1;
	line=checkF1.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))
		{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			ErrorMessage="probleme de connexion à la Database PC_InsertIntoDBLog";
			ErrorCode="1";
			exit(1);
		}
		else{
			ErrorMessage="Success connexion DB PC_InsertIntoDBLog";
			ErrorCode="0";
		}
	}
	//END INSERT INTO DB		
	
	linenum++;
	stringstream checkF2;
	checkF2<<linenum;	
	timestamp=FormatDate();
	variable="Open_FILE";
	value=value_checkFOUT2;
	type="TEXT";
	comment="ouverture fichier "+Fileout2;
	line=checkF2.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))
		{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			ErrorMessage="probleme de connexion à la Database PC_InsertIntoDBLog";
			ErrorCode="1";
			exit(1);
		}
		else{
			ErrorMessage="Success connexion DB PC_InsertIntoDBLog";
			ErrorCode="0";
		}
	}
	//END INSERT INTO DB	
	
	while(fin.good()){
		string ligne1="";
		string ligne2="";
		string ligne3="";
		string ligne4="";

		getline(fin,ligne1);
		getline(fin,ligne2);
		getline(fin,ligne3);
		getline(fin,ligne4);

		if ((!ligne1.empty())&&(!ligne2.empty())&&(!ligne3.empty())&&(!ligne4.empty())){ //au cas ou le fichier contient des retours chariots
			comptread++;
			stringstream ss1;
			stringstream ss2;
			stringstream ss3;
			stringstream ss4;
			ss1<<ligne1;
			ss2<<ligne2;
			ss3<<ligne3;
			ss4<<ligne4;

			string header1 = ss1.str();
			string seq = ss2.str();
			string header2 = ss3.str();
			string quality = ss4.str();
			if (maplines.count(seq)>0)
			{
				maplines[seq]++;
				int score =0;
				if(CasavaVersion.compare("1.7")==0)
				{
					score = membre->calcPhredScore(quality);
				}
				if(CasavaVersion.compare("1.8")==0)
				{
					score = membre->calcPhredScoreSanger(quality);
				}
				// Selection of the redundant sequence with the best phred score.
				// Add of the redundancy position of the best sequence in the hash 
				if (scorelines.count(seq)>0){
					if (score > scorelines[seq]){
						scorelines[seq] = score;
						indicelines[seq] = comptread;
					}
				}
				else{
					scorelines[seq] = score;
					indicelines[seq] = comptread;
				}
			}
			else {maplines[seq]=1;}
		}
	}
	fin.close();
	for (map<string,int>::iterator p=maplines.begin(); p!=maplines.end();++p){
		if ((*p).second>=2){
			stringstream ss;
			ss<<(*p).second;
			ss<<" ";
			ss<<(*p).first;
			string liststring = ss.str();
			listlines.push_back(liststring);
		}
	}

	listlines.sort(compare_readrep);
	for (list<string>::iterator p=listlines.begin(); p!=listlines.end();++p){
		fout2 << (*p) << endl;
	}

	fout2.close();
	maplines.clear();
	listlines.clear();
	// Creation of the index list of the sequence of interest  
	for (map<string,int>::iterator p=indicelines.begin(); p!=indicelines.end();++p){
		listindex.push_back((*p).second);
	}

	int count=0;
	//writing of uniq read files (with the best redundancy position)
	fin.open((AnalyseREP+"/"+Fileout1).c_str());
	comptread=0;
	count=0;
	while(fin.good()){
		string ligne1="";
		string ligne2="";
		string ligne3="";
		string ligne4="";
		comptread++;
		getline(fin,ligne1);
		if(!fin.eof()){
			getline(fin,ligne2);
			getline(fin,ligne3);
			getline(fin,ligne4);
			if(indicelines.count(ligne2)>0 &&(comptread==indicelines[ligne2])){
				fout<<ligne1<<endl;
				fout<<ligne2<<endl;
				fout<<ligne3<<endl;
				fout<<ligne4<<endl;
				count++;	
			}
			if(indicelines.count(ligne2)==0){
				fout<<ligne1<<endl;
				fout<<ligne2<<endl;
				fout<<ligne3<<endl;
				fout<<ligne4<<endl;
				count++;
			}
		}
	}
	cout<<"nb sequences uniques "<<Fileout1<<"\t"<<count<<endl;
	char nbUnique[20];
	sprintf (nbUnique, "%d", count);

	linenum++;
	stringstream seq1u;
	seq1u<<linenum;	
	timestamp=FormatDate();
	variable="SEQUENCE1_ANAL";
	value=Fileout1;
	type="TEXT";
	comment="nom de sequence nettoyee";
	line=seq1u.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	SampleData["SEQUENCE1_ANAL"]=Fileout1;
	
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))
		{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			ErrorMessage="probleme de connexion à la Database PC_InsertIntoDBLog";
			ErrorCode="1";
			exit(1);
		}
		else{
			ErrorMessage="Success connexion DB PC_InsertIntoDBLog";
			ErrorCode="0";
		}
	}
	//END INSERT INTO DB

	linenum++;
	stringstream nbu;
	nbu<<linenum;	
	timestamp=FormatDate();
	variable="NBunique";
	value=nbUnique;
	type="INTEGER";
	comment="nombre de sequences uniques";
	line=nbu.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))
		{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			ErrorMessage="probleme de connexion à la Database PC_InsertIntoDBLog";
			ErrorCode="1";
			exit(1);
		}
		else{
			ErrorMessage="Success connexion DB PC_InsertIntoDBLog";
			ErrorCode="0";
		}
	}
	//END INSERT INTO DB	
	
	fin.close();
	fout.close();

	//Writing in the file DBData
	fout5.open(FileDB.c_str());
	if(fout5.is_open()){
		ErrorMessage="Success";
		ErrorCode="0";
	}
	else{
		ErrorMessage="probleme à l'ouverture du fichier /inra/Jobs_server/DBData.tmp";
		ErrorCode=1;
		cout<<ErrorCode<<endl;
		cout<<ErrorMessage<<endl;
		exit(1);
	}
	for(map<string,string>::iterator si=SampleData.begin();si!= SampleData.end();si++){
		fout5<<(*si).first<<" "<<(*si).second<<" "<<endl;
	}
	fout5.close();

	// INSERTION INTO DATABASE //

	linenum++;
	stringstream endT;
	endT<<linenum;	
	string timestamp_end=FormatDate();
	cout<<"FIN DU STEP\t"<< timestamp_end<<endl;

	variable="TIME_END";
	value=timestamp_end;
	type="DATETYPE";
	comment="Temps de fin de step";
	line=endT.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";

	
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))
		{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			ErrorMessage="probleme de connexion à la Database PC_InsertIntoDBLog";
			ErrorCode="1";
			exit(1);
		}
		else{
			ErrorMessage="Success connexion DB PC_InsertIntoDBLog";
			ErrorCode="0";
		}
	}
	//END INSERT INTO DB
	mysql_close(connect);
}