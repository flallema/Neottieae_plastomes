//============================================================================
// Name        : Trimming.cpp
// Author      : Maria Tchoumakov & Mathieu Charles
// Version     : 1.7.6
// Description : Trimming of the input fastq file 
// UPDATE (by Elodie Marquand) - 20/04/2015: new version. 
// UPDATE (by Sivasangari Nandy) - 07/03/2016 : new version. Insertion in DB (var\tval)
//
// INPUT : FASTQ file(s), trimming parameter file
// OUTPUT : trimmed uniq FASTQ file(s)
//============================================================================

using namespace std;

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

void trimmingpairend(string,string,map<string,string>,map<string,string>,map<string,string>);
map<string,string> ReadDataFile(string);
string FormatDate();

class Test{
    public:
       int calcPhredScore(string);
       int calcPhredScoreSanger(string);
       
       int calcPhredScoreIndex(string,int);
       int calcPhredScoreIndexSanger(string,int);
       
       float calcPhredScoreMoyen(string,int, int);
       float calcPhredScoreMoyenSanger(string,int, int);
       
       string extractString(string,int, int);
       int num;
       string header1;
       string seq;
       string header2;
       string quality;
       int score;
       int indice;

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

int Test::calcPhredScoreSanger(string scoreline){
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
// timestamp in DB
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

map<string,string>ReadDataFile(string FileName){
	string line;
	ifstream fin;
	vector<string> words;
	map<string,string>SampleData;
	map<string,string>::iterator i;
	pair<map<string,string>::iterator,bool> ret;
	fin.open(FileName.c_str());
	if(fin.is_open()){
		cout<<"fichier "<<FileName<<" ouvert\n";
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

int main(int argc,char*argv[]){

	string line;
	int count=0;
	string Id_Processus=argv[1];
	string FileName=argv[2];
	cout<<"FileName\t"<<FileName<<endl;
	string ParamFile=argv[3];
	cout<<"ParamFile\t" <<ParamFile<<endl;
	string FileConnectTmp=argv[4];
	cout<<"FileConnectTmp\t"<<FileConnectTmp<<endl;
	
	map<string,string>Connector=ReadDataFile(FileConnectTmp);
	map<string,string>parametres=ReadDataFile(ParamFile);
	map<string,string>Sample_Data=ReadDataFile(FileName);

    trimmingpairend(Id_Processus,FileName,parametres,Sample_Data,Connector);
	return 0;
}

void trimmingpairend(string Id_Processus,string FileName,map<string,string>parametres,map<string,string>Sample_Data,map<string,string>Connector){

	int linenum=0;
	string step="";
	string timestamp="";
	string variable="";
	string value="";
	string type="";
	string comment="";
	string line="";
	string step_number="";
	string query="";	
	string ErrorCode="";
	string ErrorMessage="";
	MYSQL *connect;
	MYSQL *status;
	MYSQL_ROW row;
	MYSQL_RES *result;	
	string server=Connector["host"];
	string database=Connector["dbname"];
	string user=Connector["user"];
	string pass=Connector["password"];
	string timestamp_deb=FormatDate();
	cout<<"DEBUT DU STEP\t"<< timestamp_deb<<endl;
	connect=mysql_init(NULL);
	char* cmd;
	string queryStep="call PC_GetNextStep(@idProcessus:='"+Id_Processus+"');";
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
	
	ifstream fin1;
	ifstream fin2;
	ofstream fout1;
	ofstream fout2;
	ofstream fileEliminatedRead1;
	ofstream fileEliminatedRead2;	
	ofstream fileEliminatedRead1afterN;
	ofstream fileEliminatedRead2afterN;		
	
	ofstream fileStats;	
	unsigned int comptread=0;
	unsigned int comptjoin=0;
	unsigned int compttrimmed=0;
	map<int,int>comptN;
	int comptN1=0;
	int comptN2=0;
	Test* membre = new Test;
	int nbN=0;
	unsigned int  nbRead1elimLmin=0;
	unsigned int  nbRead2elimLmin=0;
	unsigned int  nbRead1elimAfterN=0;
	unsigned int  nbRead2elimAfterN=0;
	
	string RepQC=Sample_Data["REPPARENT"];
	string RepDecomp=Sample_Data["REPDECOMP"];
	string AnalyseREP=Sample_Data["REPANAL"];
	string Seq1=Sample_Data["SEQUENCE1_INIT"];
	string Seq2=Sample_Data["SEQUENCE2_INIT"];
	string CasavaVersion=Sample_Data["Casava_Version"];
	
	string Fileout1="trimmed_"+Seq1;
	string Fileout2="trimmed_"+Seq2;
	string fileEliminatedRead1name="eliminatedRead1_"+Seq1;
	string fileEliminatedRead2name="eliminatedRead2_"+Seq2;
	string fileEliminatedRead1afterNname="eliminatedRead1AfterN_"+Seq1;
	string fileEliminatedRead2afterNname="eliminatedRead2AfterN_"+Seq2;	
	string fileStatsName="stats_"+Seq1;
		
	// Get the Trimming parameter from input file
	int cyclemax=atoi(parametres["taille_sequences"].c_str());
	int Qminibase=atoi(parametres["qualite_mini"].c_str());
	int Qminimoy=atoi(parametres["moyenne_mini"].c_str());
	int nbbasemin=atoi(parametres["taille_mini"].c_str());
	int Qmini=atoi(parametres["limit_score(>=)"].c_str());
	int filtreN=atoi(parametres["max_nbN_accepte"].c_str());
	
	if(cyclemax ==0 || Qminibase==0 || Qminimoy==0 || nbbasemin==0 || Qmini==0)
	{
		ErrorCode="1";
		ErrorMessage="Le contenu/syntaxe du fichier de parametre n'est pas correct";
	}
	
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream si;
	si<<linenum;
	cout<<step_number<<endl;
	mysql_free_result(result);
	step="Trimming";
	timestamp=FormatDate();
	variable="TIME_START";
	value=timestamp_deb;
	type="DATETYPE";
	comment="Temps du debut du step";
	line=si.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
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
	
	cout<<"cyclemax\t"<<cyclemax<<endl;
	cout<<"qualite min par base\t"<<Qminibase<<endl;	
	cout<<"qualite moy min sequence\t"<<Qminimoy<<endl;	
	cout<<"nbbasemin\t"<<nbbasemin<<endl;	
	cout<<"qualite limite score>=\t"<<Qmini<<endl;	
	cout<<"N accepte\t"<<filtreN<<endl;	

	string value_def1="";
	string value_def2="";
	fout1.open((AnalyseREP+'/'+Fileout1).c_str());
	if(fout1.is_open()){
		cout<<"Ouverture fichier trimmed_R1 OK\t"<<Fileout1<<endl;	
		value_def1="Ouverture fichier trimmed_R1 OK : "+Fileout1+"\n";
	}
	else{
		cout<<"ERROR : probleme d'ouverture "<<Fileout1<<endl;
		value_def1="ERROR : probleme d'ouverture "+Fileout1+"\n";
	}
	fout2.open((AnalyseREP+'/'+Fileout2).c_str());
	if(fout2.is_open()){
		cout<<"Ouverture fichier trimmed_R2 OK\t"<<Fileout2<<endl;
		value_def2="Ouverture fichier trimmed_R2 OK : "+Fileout2+"\n";
	}
	else{
		cout<<"ERROR : probleme d'ouverture "<<Fileout2<<endl;
		value_def2="ERROR : probleme d'ouverture "+Fileout2+"\n";
	}
	
	// Check trimmed uniq files opening
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream sif;
	sif<<linenum;	
	timestamp=FormatDate();
	variable="CHECK_FILE";
	value=value_def1;
	type="TEXT";
	comment="controle ouverture de fichier";
	line=sif.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
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
	
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream sif2;
	sif2<<linenum;	
	timestamp=FormatDate();
	variable="CHECK_FILE";
	value=value_def2;
	type="TEXT";
	comment="controle ouverture de fichier";
	line=sif2.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
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

	// Read == uniq sequences
	string value_defU1="";
	string value_defU2="";	
	string Read1=RepDecomp+'/'+Seq1;
	string Read2=RepDecomp+'/'+Seq2;
	fin1.open(Read1.c_str());
	if(fin1.is_open()){
		cout<<"Ouverture fichier uniq_R1 OK\t"<<Seq1<<endl;
		value_defU1="Ouverture fichier uniq_R1 OK : "+Seq1+"\n";		
	}
	else{
		cout<<"ERROR : probleme d'ouverture "<<Seq1<<endl;
		value_defU1="ERROR : probleme d\'ouverture "+Seq1+"\n";
	}
	fin2.open(Read2.c_str());
	if(fin2.is_open()){
		cout<<"Ouverture fichier uniq_R2 OK\t"<<Seq2<<endl;
		value_defU2="Ouverture fichier uniq_R2 OK : "+Seq2+"\n";
	}
	else{
		cout<<"ERROR : probleme d'ouverture "<<Seq2<<endl;
		value_defU2="ERROR : probleme d\'ouverture "+Seq2+"\n";
	}
	
	// Check of uniq read files opening
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream sifU1;
	sifU1<<linenum;	
	timestamp=FormatDate();
	variable="CHECK_FILE";
	value=value_defU1;
	type="TEXT";
	comment="controle ouverture de fichier";
	line=sifU1.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
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
	
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream sifU2;
	sifU2<<linenum;	
	timestamp=FormatDate();
	variable="CHECK_FILE";
	value=value_defU2;
	type="TEXT";
	comment="controle ouverture de fichier";
	line=sifU2.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
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

	string value_elimR1="";
	string value_elimR2="";
	string value_elimR1_afterN="";
	string value_elimR2_afterN="";
	string value_statU="";	
	
	fileEliminatedRead1.open((AnalyseREP+'/'+fileEliminatedRead1name).c_str());
	if(fileEliminatedRead1.is_open()){
		cout<<"Ouverture fichier Eliminated_R1 OK\t"<<fileEliminatedRead1name<<endl;	
		value_elimR1="Ouverture fichier Eliminated_R1 OK : "+fileEliminatedRead1name+"\n";
	}
	else{
		cout<<"ERROR : probleme d'ouverture  "<<fileEliminatedRead1name<<endl;
		value_elimR1="ERROR : probleme d'ouverture  "+fileEliminatedRead1name+"\n";	
	}
	fileEliminatedRead2.open((AnalyseREP+'/'+fileEliminatedRead2name).c_str());
	if(fileEliminatedRead2.is_open()){
		cout<<"Ouverture fichier Eliminated_R2 OK\t"<<fileEliminatedRead2name<<endl;	
		value_elimR2="Ouverture fichier Eliminated_R2 OK : "+fileEliminatedRead2name+"\n";
	}
	else{
		cout<<"ERROR : probleme d'ouverture  "<<fileEliminatedRead2name<<endl;
		value_elimR2="ERROR : probleme d'ouverture  "+fileEliminatedRead2name+"\n";			
	}
	fileEliminatedRead1afterN.open((AnalyseREP+'/'+fileEliminatedRead1afterNname).c_str());
	if(fileEliminatedRead1afterN.is_open()){
		cout<<"Ouverture fichier eliminatedR1AfterN OK\t"<<fileEliminatedRead1afterNname<<endl;	
		value_elimR1_afterN="Ouverture fichier eliminatedR1AfterN OK : "+fileEliminatedRead1afterNname+"\n";
	}
	else{
		cout<<"ERROR : probleme d'ouverture  "<<fileEliminatedRead1afterNname<<endl;
		value_elimR1_afterN="ERROR : probleme d'ouverture  "+fileEliminatedRead1afterNname+"\n";					
	}
	fileEliminatedRead2afterN.open((AnalyseREP+'/'+fileEliminatedRead2afterNname).c_str());
	if(fileEliminatedRead2afterN.is_open()){
		cout<<"Ouverture fichier eliminatedR2AfterN OK\t"<<fileEliminatedRead2afterNname<<endl;	
		value_elimR2_afterN="Ouverture fichier eliminatedR2AfterN OK : "+fileEliminatedRead2afterNname+"\n";		
	}
	else{
		cout<<"ERROR : probleme d'ouverture  "<<fileEliminatedRead2afterNname<<endl;
		value_elimR2_afterN="ERROR : probleme d'ouverture  "+fileEliminatedRead2afterNname+"\n";							
	}

	fileStats.open((AnalyseREP+'/'+fileStatsName).c_str());
	if(fileStats.is_open()){
		cout<<"Ouverture fichier stat OK\t"<<fileStatsName<<endl;	
		value_statU="Ouverture fichier stat OK : "+fileStatsName+"\n";				
	}
	else{
		cout<<"ERROR : probleme d'ouverture  "<<fileStatsName<<endl;
		value_statU="ERROR : probleme d'ouverture  "+fileStatsName+"\n";									
	}
	
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream ver1;
	ver1<<linenum;	
	timestamp=FormatDate();
	variable="CHECK_FILE";
	value=value_elimR1;
	type="TEXT";
	comment="controle ouverture de fichier";
	line=ver1.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
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
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream ver2;
	ver2<<linenum;	
	timestamp=FormatDate();
	variable="CHECK_FILE";
	value=value_elimR2;
	type="TEXT";
	comment="controle ouverture de fichier";
	line=ver2.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
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
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream ver1N;
	ver1N<<linenum;	
	timestamp=FormatDate();
	variable="CHECK_FILE";
	value=value_elimR1_afterN;
	type="TEXT";
	comment="controle ouverture de fichier";
	line=ver1N.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
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
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream ver2N;
	ver2N<<linenum;	
	timestamp=FormatDate();
	variable="CHECK_FILE";
	value=value_elimR2_afterN;
	type="TEXT";
	comment="controle ouverture de fichier";
	line=ver2N.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
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
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream statU;
	statU<<linenum;	
	timestamp=FormatDate();
	variable="CHECK_FILE";
	value=value_statU;
	type="TEXT";
	comment="controle ouverture de fichier";
	line=statU.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
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
	
	int countlines=0;
	fin1.clear();
	fin1.seekg(0);
	fin2.clear();
	fin2.seekg(0);
	while((fin1.good())&&((fin2.good()))){
 		if(!fin1.eof() && !fin2.eof()){
			countlines++;
		} 
		comptread++;
		string r1line1;
		string r1line2;
		string r1line3;
		string r1line4;

		string r2line1;
		string r2line2;
		string r2line3;
		string r2line4;

		getline(fin1,r1line1);
		getline(fin1,r1line2);
		getline(fin1,r1line3);
		getline(fin1,r1line4);

		getline(fin2,r2line1);
		getline(fin2,r2line2);
		getline(fin2,r2line3);
		getline(fin2,r2line4);
		//SN add: 
		int maxlen=0;

		// SN: Particular case of sequence headers 
		// @HWI-ST1194:55:C11P5ACXX:4:1101:1219:2151 2:N:0:
		// @HISEQ8:C3FVWACXX:6:1101:1489:1921/1

		if ((r1line1.find('/')) < 80 )
		{
			maxlen=r1line1.find('/');
		}
		if ((r1line1.find(' ')) < 80)
		{
			maxlen=r1line1.find(' ');
		}
		if ((!r1line1.empty())&&(r1line1.substr(0,r1line1.size()-1) == r2line1.substr(0,r2line1.size()-1))){

			stringstream r1ss1;
			stringstream r1ss2;
			stringstream r1ss3;
			stringstream r1ss4;

			stringstream r2ss1;
			stringstream r2ss2;
			stringstream r2ss3;
			stringstream r2ss4;

			r1ss1<<r1line1;
			r1ss2<<r1line2;
			r1ss3<<r1line3;
			r1ss4<<r1line4;

			r2ss1<<r2line1;
			r2ss2<<r2line2;
			r2ss3<<r2line3;
			r2ss4<<r2line4;

			string r1header1 = r1ss1.str();
			string r1seq = r1ss2.str();
			string r1header2 = r1ss3.str();
			string r1quality = r1ss4.str();

			string r2header1 = r2ss1.str();
			string r2seq = r2ss2.str();
			string r2header2 = r2ss3.str();
			string r2quality = r2ss4.str();
			
			int z1=0;
			for (z1=cyclemax-1;z1>=0;--z1)
			{
				int r1qscorebase;
				float r1qscoremoy;
				if(CasavaVersion.compare("1.8")==0)	
				{
					r1qscorebase = membre->calcPhredScoreIndexSanger(r1quality,z1);
					r1qscoremoy = membre->calcPhredScoreMoyenSanger(r1quality, 0, z1);
				}
				else
				{
					r1qscorebase = membre->calcPhredScoreIndex(r1quality,z1);
					r1qscoremoy = membre->calcPhredScoreMoyen(r1quality, 0, z1);
				}
				
			
				if (r1qscorebase > Qmini){
					if (r1qscorebase>=Qminibase){
						break;
					}
					else {
						if (r1qscoremoy>=Qminimoy){
							break;
						}
						else{
						}
					}
				}
				
				
			}
			
			int z2=0;
			for (z2=cyclemax-1;z2>=0;--z2)
			{
				int r2qscorebase;
				float r2qscoremoy;
				if(CasavaVersion.compare("1.8")==0)	
				{
					r2qscorebase = membre->calcPhredScoreIndexSanger(r2quality,z2);
					r2qscoremoy = membre->calcPhredScoreMoyenSanger(r2quality, 0, z2);
					
				}
				else
				{	
					r2qscorebase = membre->calcPhredScoreIndex(r2quality,z2);
					r2qscoremoy = membre->calcPhredScoreMoyen(r2quality, 0, z2);
				}
							
				if (r2qscorebase > Qmini){
					if (r2qscorebase>=Qminibase){
						break;
					}
					else {
						if (r2qscoremoy>=Qminimoy){
							break;
						}
						else{
						}
					}
				}
			}
	
		
			bool write = true; 
			if (z1>=nbbasemin-1){
				if (z2>=nbbasemin-1){
					nbN=0;
					string read1trimme=membre->extractString(r1seq,0,z1);
					string read2trimme=membre->extractString(r2seq,0,z2);
					for(string::iterator si=read1trimme.begin();si!=read1trimme.end();si++){
						if(*si=='N'){
							nbN++;
						}
						if(nbN>filtreN){

							write=false;
							comptN[comptread]=nbN;
							comptN1++;
						}
					}
					
					if(write==false){
						//read1 eliminated
						nbRead1elimAfterN++;
						fileEliminatedRead1afterN << r1header1<<endl; //ELM addition
						fileEliminatedRead1afterN << membre->extractString(r1seq,0,z1) <<endl;//ELM addition
						fileEliminatedRead1afterN << r1header2 <<endl;//ELM addition
						fileEliminatedRead1afterN << membre->extractString(r1quality,0,z1) <<endl;//ELM addition	
					}
					if(write==true){
						nbN=0;
						for(string::iterator si=read2trimme.begin();si!=read2trimme.end();si++){
							if(*si=='N'){
								nbN++;
								if(nbN>filtreN){
									write=false;
									comptN[comptread]=nbN;
									comptN2++;
								}
							}
						}
						if(write==false){
							//read2 eliminated
							nbRead2elimAfterN++;
							fileEliminatedRead2afterN << r2header1<<endl; //ELM addition
							fileEliminatedRead2afterN << membre->extractString(r2seq,0,z2) <<endl;//ELM addition
							fileEliminatedRead2afterN << r2header2 <<endl;//ELM addition
							fileEliminatedRead2afterN << membre->extractString(r2quality,0,z2) <<endl;//ELM addition	
						}
					}
					if (write==true){
						fout1 << r1header1<<endl;
						fout1 << membre->extractString(r1seq,0,z1) <<endl;
						fout1 << r1header2 <<endl;
						fout1 << membre->extractString(r1quality,0,z1) <<endl;
						fout2<< r2header1 << endl;
						fout2<< membre->extractString(r2seq,0,z2) <<endl;
						fout2<< r2header2 << endl;
						fout2<< membre->extractString(r2quality,0,z2) <<endl;
						comptjoin++;
					}
				}else{
					nbRead2elimLmin++;
					fileEliminatedRead2 << r2header1<<endl; //ELM addition
					fileEliminatedRead2 << membre->extractString(r2seq,0,z2) <<endl;//ELM addition
					fileEliminatedRead2 << r2header2 <<endl;//ELM addition
					fileEliminatedRead2 << membre->extractString(r2quality,0,z2) <<endl;//ELM addition	
					compttrimmed++;
				}
			}
			else {
				nbRead1elimLmin++;
				fileEliminatedRead1 << r1header1<<endl; //ELM addition
				fileEliminatedRead1 << membre->extractString(r1seq,0,z1) <<endl;//ELM addition
				fileEliminatedRead1 << r1header2 <<endl;//ELM addition
				fileEliminatedRead1 << membre->extractString(r1quality,0,z1) <<endl;//ELM addition	
				compttrimmed++;
				
			}
		}
		else {
		}
	}
	countlines=countlines-1;
	cout << "PE nb sequences uniques"<<"\t"<<countlines<<endl;

	
	//29/04/2014 : check of files closing 
	string value_r1_close="";
	string value_r2_close="";
	string value_elimR1_close="";
	string value_elimR2_close="";
	string file1_close="";
	string file2_close="";
	string value_elimR1_afterN_close="";
	string value_elimR2_afterN_close="";
	string value_statU_close="";
		
	fin1.close();
	if(fin1.is_open()){
		cout<<"ERROR : probleme de fermeture  "<<Read1<<endl;
		value_r1_close="ERROR : probleme de fermeture "+Read1+"\n";	
	}else{
		cout<<"Fermeture fichier OK\t"<<Read1<<endl;	
		value_r1_close="Fermeture fichier OK : "+Read1+"\n";			
	}
	fin2.close();
	if(fin2.is_open()){
		cout<<"ERROR : probleme de fermeture  "<<Read2<<endl;
		value_r2_close="ERROR : probleme de fermeture "+Read2+"\n";			
	}else{
		cout<<"Fermeture fichier OK\t"<<Read2<<endl;
		value_r2_close="Fermeture fichier OK : "+Read2+"\n";					
	}
	fout1.close();
	if(fout1.is_open()){
		cout<<"ERROR : probleme de fermeture  "<<Fileout1<<endl;
		file1_close="ERROR : probleme de fermeture "+Fileout1+"\n";							
		
	}else{
		cout<<"Fermeture fichier OK\t"<<Fileout1<<endl;	
		file1_close="Fermeture fichier OK : "+Fileout1+"\n";									
	}
	fout2.close();
	if(fout2.is_open()){
		cout<<"ERROR : probleme de fermeture  "<<Fileout2<<endl;
		file2_close="ERROR : probleme de fermeture "+Fileout2+"\n";											
	}else{
		cout<<"Fermeture fichier OK\t"<<Fileout2<<endl;
		file2_close="Fermeture fichier OK : "+Fileout2+"\n";						
	}
	
	fileEliminatedRead1.close();
	if(fileEliminatedRead1.is_open()){
		cout<<"ERROR : probleme de fermeture  "<<fileEliminatedRead1name<<endl;
		value_elimR1_close="ERROR : probleme de fermeture "+fileEliminatedRead1name+"\n";													
	}else{
		cout<<"Fermeture fichier OK\t"<<fileEliminatedRead1name<<endl;	
		value_elimR1_close="Fermeture fichier OK : "+fileEliminatedRead1name+"\n";															
	}
	fileEliminatedRead2.close();
	if(fileEliminatedRead2.is_open()){
		cout<<"ERROR : probleme de fermeture  "<<fileEliminatedRead2name<<endl;
		value_elimR2_close="ERROR : probleme de fermeture "+fileEliminatedRead2name+"\n";															
	}else{
		cout<<"Fermeture fichier OK\t"<<fileEliminatedRead2name<<endl;	
		value_elimR2_close="Fermeture fichier OK : "+fileEliminatedRead2name+"\n";																	
	}
	fileEliminatedRead1afterN.close();
	if(fileEliminatedRead1afterN.is_open()){
		cout<<"ERROR : probleme de fermeture  "<<fileEliminatedRead1afterNname<<endl;
		value_elimR1_afterN_close="ERROR : probleme de fermeture "+fileEliminatedRead1afterNname+"\n";																	
	}else{
		cout<<"Fermeture fichier OK\t"<<fileEliminatedRead1afterNname<<endl;	
		value_elimR1_afterN_close="Fermeture fichier OK : "+fileEliminatedRead1afterNname+"\n";																			
	}
	fileEliminatedRead2afterN.close();
	if(fileEliminatedRead2afterN.is_open()){
		cout<<"ERROR : probleme de fermeture  "<<fileEliminatedRead2afterNname<<endl;
		value_elimR2_afterN_close="ERROR : probleme de fermeture "+fileEliminatedRead2afterNname+"\n";																			
	}else{
		cout<<"Fermeture fichier OK\t"<<fileEliminatedRead2afterNname<<endl;	
		value_elimR2_afterN_close="Fermeture fichier OK : "+fileEliminatedRead2afterNname+"\n";																					
	}
	
	fileStats << "Seq _ : " << Seq1 << endl;	
	fileStats << "nb read 1 elimines (Longueur mini) : " << nbRead1elimLmin << endl;	
	fileStats << "nb read 2 elimines (Longueur mini) : " << nbRead2elimLmin << endl;	
	fileStats << "nb read 1 elimines (Nb N) : " << nbRead1elimAfterN << endl;	
	fileStats << "nb read 2 elimines (Nb N) : " << nbRead2elimAfterN << endl;	
	fileStats << "Couples de reads PE elimines par trimming (seuil de taille mini) : " << compttrimmed<< endl;
	fileStats << "Couples de reads PE elimines  > Nmax  : " << comptN.size() << endl;
	fileStats << "PE nb sequences apres nettoyage : " << comptjoin << endl;
	
	fileStats.close();
	if(fileStats.is_open()){
		cout<<"ERROR : probleme de fermeture  "<<fileStatsName<<endl;
		value_statU_close="ERROR : probleme de fermeture "+fileStatsName+"\n";																					
	}else{
		cout<<"Fermeture fichier OK\t"<<fileStatsName<<endl;	
		value_statU_close="Fermeture fichier OK : "+fileStatsName+"\n";																					
	}
	
	cout << "nb read 1 elimines (Longueur mini)\t" << nbRead1elimLmin << endl;	
	cout << "nb read 2 elimines (Longueur mini)\t" << nbRead2elimLmin << endl;	
	cout << "nb read 1 elimines (Nb N)\t" << nbRead1elimAfterN << endl;	
	cout << "nb read 2 elimines (Nb N)\t" << nbRead2elimAfterN << endl;	
	cout << "Couples de reads PE elimines par trimming (seuil de taille mini)\t" << compttrimmed<< endl;
	cout << "Couples de reads PE elimines  > Nmax\t" << comptN.size() << endl;
	cout << "PE nb sequences apres nettoyage\t" << comptjoin << endl;

	// INSERTION INTO DATABASE //
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
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
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
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream r1close;
	r1close<<linenum;	
	timestamp=FormatDate();
	variable="CHECK_FILE";
	value=value_r1_close;
	type="TEXT";
	comment="controle fermeture de fichier";
	line=r1close.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
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
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream r2close;
	r2close<<linenum;	
	timestamp=FormatDate();
	variable="CHECK_FILE";
	value=value_r2_close;
	type="TEXT";
	comment="controle fermeture de fichier";
	line=r2close.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
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
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream r1elmclose;
	r1elmclose<<linenum;	
	timestamp=FormatDate();
	variable="CHECK_FILE";
	value=value_elimR1_close;
	type="TEXT";
	comment="controle fermeture de fichier";
	line=r1elmclose.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
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
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream r2elmclose;
	r2elmclose<<linenum;	
	timestamp=FormatDate();
	variable="CHECK_FILE";
	value=value_elimR2_close;
	type="TEXT";
	comment="controle fermeture de fichier";
	line=r2elmclose.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
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
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream f1close;
	f1close<<linenum;	
	timestamp=FormatDate();
	variable="CHECK_FILE";
	value=file1_close;
	type="TEXT";
	comment="controle fermeture de fichier";
	line=f1close.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
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
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream f2close;
	f2close<<linenum;	
	timestamp=FormatDate();
	variable="CHECK_FILE";
	value=file2_close;
	type="TEXT";
	comment="controle fermeture de fichier";
	line=f2close.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
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
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream elmNclose;
	elmNclose<<linenum;	
	timestamp=FormatDate();
	variable="CHECK_FILE";
	value=value_elimR1_afterN_close;
	type="TEXT";
	comment="controle fermeture de fichier";
	line=elmNclose.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
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
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream elm2Nclose;
	elm2Nclose<<linenum;	
	timestamp=FormatDate();
	variable="CHECK_FILE";
	value=value_elimR2_afterN_close;
	type="TEXT";
	comment="controle fermeture de fichier";
	line=elm2Nclose.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
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
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream statUclose;
	statUclose<<linenum;	
	timestamp=FormatDate();
	variable="CHECK_FILE";
	value=value_statU_close;
	type="TEXT";
	comment="controle fermeture de fichier";
	line=statUclose.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
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
	
	char Tcount1[20];
	sprintf (Tcount1, "%d", nbRead1elimLmin);	
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream nb1;
	nb1<<linenum;	
	timestamp=FormatDate();
	variable="NB_read1_elimine_longueurMIN";
	value=Tcount1;
	type="INTEGER";
	comment="nb read 1 elimines (Longueur mini)";
	line=nb1.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
		
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
	
	char Tcount[20];
	sprintf (Tcount, "%d", nbRead2elimLmin);	
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream nb2;
	nb2<<linenum;	
	timestamp=FormatDate();
	variable="NB_read2_elimine_longueurMIN";
	value=Tcount;
	type="INTEGER";
	comment="nb read 2 elimines (Longueur mini)";
	line=nb2.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
		
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
	
	char Tcountv[20];
	sprintf (Tcountv, "%d", nbRead1elimAfterN);		
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream nbn1;
	nbn1<<linenum;	
	timestamp=FormatDate();
	variable="NB_read1_elimine_Nb_N";
	value=Tcountv;
	type="INTEGER";
	comment="nb read 1 elimines (Nb N)";
	line=nbn1.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
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
	
	char Tcountv2[20];
	sprintf (Tcountv2, "%d", nbRead2elimAfterN);		
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream nbn2;
	nbn2<<linenum;	
	timestamp=FormatDate();
	variable="NB_read2_elimine_Nb_N";
	value=Tcountv2;
	type="INTEGER";
	comment="nb read 2 elimines (Nb N)";
	line=nbn2.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
		
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
		
	char Tcounttrim[20];
	sprintf (Tcounttrim, "%d", compttrimmed);	
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream pee;
	pee<<linenum;	
	timestamp=FormatDate();
	variable="NB_readsPE_elimine_par_Trimming";
	value=Tcounttrim;
	type="INTEGER";
	comment="reads PE elimines par trimming (seuil de taille mini)";
	line=pee.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	
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

	char TcountN[20];
	sprintf (TcountN, "%d", comptN.size());		
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream pec;
	pec<<linenum;	
	timestamp=FormatDate();
	variable="NB_readsPE_elimine";
	value=TcountN;
	type="INTEGER";
	comment="Couples de reads PE elimines  > Nmax";
	line=pec.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
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

	char Tcountjoin[20];
	sprintf (Tcountjoin, "%d", comptjoin);			
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream nsn;
	nsn<<linenum;	
	timestamp=FormatDate();
	variable="NB_seq_apres_nettoyage";
	value=Tcountjoin;
	type="INTEGER";
	comment="PE nb sequences apres nettoyage";
	line=nsn.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
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
    query="call PC_InsertIntoDBLog('"+Id_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
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


