//============================================================================
// Name        : QualityControl.cpp
// Author      : Maria Tchoumakov & Mathieu Charles
// Version     : 1.7.6
// Description : Quality control of the input fastq file 
// Amount of ATGCN; Data table creation for boxplot
//
// UPDATE (by Sivasangari Nandy) - 07/03/2016 : new version. Insertion in DB (var\tval)
//
// INPUT : FASTQ file
// OUTPUT : files *ATGCN_stats.txt & *boxplotvalues.txt
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
#include <string.h>
#include <sstream>
#include <algorithm>
#include <map>
#include <list>
#include <exception>
#include <mysql.h>
#include <time.h>

void createStat (string,string,map<string,string>,map<string,string>); 

int ControlRead(string,string,map<string,string>,map<string,string>);
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
       
       Test ()
       {
       //    cout << "pouet";
       }
       
       Test (string laseq)
       {
           seq=laseq;
       }
       
       Test (int Num, string Header1, string Seq, string Header2, string Quality, int Indice)
       {
           num = Num;
           header1 = Header1;
           seq=Seq;
           header2 = Header2;
           quality = Quality;
           score=this->calcPhredScore(Quality);
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

int Test::calcPhredScoreSanger(string scoreline){ //calcul score PHRED SANGER (illumina 1.8 et suivant)
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

	int size = seq.size();
	char * scoretab = new char[ size ];
	strncpy( scoretab, seq.c_str(), size);
	int score = (int) scoretab[index]-64;

	delete [] scoretab;

	return score;
}

int Test::calcPhredScoreIndexSanger (string seq, int index){ //calcul score PHRED SANGER (illumina 1.8 et suivant)

	int size = seq.size();
	char * scoretab = new char[ size ];
	strncpy( scoretab, seq.c_str(), size);
	int score = (int) scoretab[index]-33;
	delete [] scoretab;
	return score;
}

int main(int argc,char*argv[]){
	string Read=argv[1];
	string Id_Processus=argv[2];
	string FileDBName=argv[3];
	string FileConnectTmp=argv[4];
	map<string,string> SampleData=ReadDataFile(FileDBName);
	map<string,string>Connector=ReadDataFile(FileConnectTmp);

	ControlRead(Read, Id_Processus,SampleData,Connector);

	createStat (Read, Id_Processus,SampleData,Connector );
	return 0;
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

	}
	else{
		cout<<1<<endl;
		cout<<"probleme à l'ouverture du fichier "<<FileName<<endl;
		exit(1);
	}
	int count=0;
	while(!fin.eof()){
		getline(fin,line);
		//cout<<line<<endl;
		count++;
		Split(words,line," ");
		SampleData.insert( pair<string,string>(words[0],words[1]) );
	}
	fin.close();
	return SampleData;
}

int ControlRead(string Read, string ID_Processus,map<string,string>SampleData,map<string,string>Connector)
{
	int control;
    int count=0;
	char seq[200];
	string line1;
	string line2;
	string line3;
	string line4;
	ifstream fin;
	ofstream fout;
	string step_name="QualityControl_DB ";
	string message="";
	string FileName="";
	string read1or2="";	
	
	string CasavaVersion="";
	string AnalyseREP="";
	string compress="";
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
	string repdecomp="";
	string repbrut="";	
	
	string ErrorCode="";
	string ErrorMessage="";
	MYSQL *connect;
	MYSQL *status;
	MYSQL_ROW row;
	MYSQL_RES *result;
	string readname="";
	string index;
	size_t found;
	size_t found2;
	string server=Connector["host"];
	string database=Connector["dbname"];
	string user=Connector["user"];
	string pass=Connector["password"];
	string timestamp_deb=FormatDate();
	connect=mysql_init(NULL);
	char* cmd;
	string queryStep="call PC_GetNextStep(@idProcessus:='"+ID_Processus+"');";
	AnalyseREP=SampleData["REPANAL"];
	
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,queryStep.c_str())){
			fprintf(stderr, "%s\n", mysql_error(connect));
			cout<<cmd<<endl;
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
		ErrorMessage=step_name+"\t probleme de connexion à la Database PC_GetNextStep";
		ErrorCode="1";
		// cout<<ErrorCode<<endl;
		cout<<ErrorMessage<<endl;
		exit(1);
	}
	cout<<"step num\t"<< step_number<<endl;
	cout<<"DEBUT DU STEP\t"<< timestamp_deb<<endl;

	// INSERTION INTO DATABASE //
	linenum++;
	stringstream si;
	si<<linenum;

	mysql_free_result(result);
	step="Quality_Control";
	timestamp=FormatDate();
	variable=Read+"_TIME_START";
	value=timestamp_deb;
	type="DATETYPE";
	comment="Temps du debut du step";
	line=si.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))	{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			exit(1);
		}
		else{

		}
	}
	//END INSERT INTO DB
	
	map<string,string>::iterator i ;	
	for(i=SampleData.begin();i!=SampleData.end();i++){
		if((*i).first.compare("REPDECOMP")==0){
			repdecomp=(*i).second;
		}
		if((*i).first.compare("REPBRUT")==0){
			repbrut=(*i).second;
		}
		
	}
	if(repdecomp.compare(repbrut)==0 ){
		if(Read.compare("Read1")==0){
				FileName=SampleData["Anal_Read1"];
				read1or2="R1"; 				
		}
		if( Read.compare("Read2")==0){
				FileName=SampleData["Anal_Read2"];
				read1or2="R2"; 
		}	
	}
	else{
		if(Read.compare("Read1")==0){
			FileName=SampleData["R1_Filename"];
			read1or2="R1"; 
		}
		if(Read.compare("Read2")==0){
			FileName=SampleData["R2_Filename"];
			read1or2="R2"; 						
		}
		CasavaVersion=SampleData["Casava_Version"];
	}
	
	string value_filename="";
	string value_decomp="";
	string value_cheminR="";
	
	cout<<"Nom du fichier de Read\t"<< FileName<<endl;
	value_filename="Nom du fichier de Read\t"+FileName+"\n";
	cout<<" Repertoire de decompression\t"<<repdecomp<<endl;
	value_decomp=" Repertoire de decompression "+repdecomp+"\n";
	string Path=repdecomp+"/"+FileName;
	cout<<"Chemin du Read\t"<<Path<<endl;
	value_cheminR=Path+"\n";
	
	//Sequence fastq file opening	

	fin.open(Path.c_str());
	string value_def="";
	if(fin.is_open())
	{
		cout<<" Fichier ControlRead "<<FileName<<" ouvert\n";
		value_def="Fichier ControlRead "+ FileName +" ouvert";
	}
	else
	{
		cout<<"probleme à l'ouverture du fichier "<<Path<<"\n";
		value_def= "probleme à l'ouverture du fichier "+Path;
	}

	// INSERTION INTO DATABASE //
	linenum++;
	stringstream nfn;
	nfn<<linenum;	
	timestamp=FormatDate();
	variable="Read_filename";
	value=value_filename;
	type="TEXT";
	comment="Nom du fichier de Read";
	line=nfn.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))
		{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			exit(1);
		}
		else{

		}
	}
	//END INSERT INTO DB
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream nvd;
	nvd<<linenum;	
	timestamp=FormatDate();
	variable="Decompression_directory";
	value=value_decomp;
	type="TEXT";
	comment="Repertoire de decompression";
	line=nvd.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))
		{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			exit(1);
		}
		else{

		}
	}
	//END INSERT INTO DB
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream nrp;
	nrp<<linenum;	
	timestamp=FormatDate();
	variable="Read_path";
	value=value_cheminR;
	type="TEXT";
	comment="Chemin du Read";
	line=nrp.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	//INSERT INTO DB
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))
		{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			exit(1);
		}
		else{

		}
	}
	//END INSERT INTO DB
	
	
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream sif;
	sif<<linenum;	
	timestamp=FormatDate();
	variable=Read+"_FILE";
	value=value_def;
	type="TEXT";
	comment="controle ouverture de fichier";
	line=sif.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))
		{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			exit(1);
		}
		else{

		}
	}
	//END INSERT INTO DB

	string value_checkF="";
	string value_checkC="";
	
	while(fin.good())
	{
		getline(fin,line1);
		if(!fin.eof())
		{
			getline(fin,line2);
			getline(fin,line3);
			getline(fin,line4);
			strcpy(seq,line2.c_str());
			count++;
			char str_count [20];
			sprintf (str_count, "%d", count);
			
			control=1;
			if((line1.empty() &&line1.compare("\n")==0)||((line2.empty() && line2.compare("\n")==0))||(line3.empty() && line3.compare("\n")==0)||(line4.empty() && line4.compare("\n")==0))
			{
				control=1;
				cout<<"lignes vides dans le fichier "<<FileName<<" position "<<count<<"\n";
				value_checkF="lignes vides dans le fichier "+FileName+" position "+str_count;
				exit;
			}
			if(line1.compare("@")==0)
			{
				control=1;
				cout<<"Erreur dans le format du nom de read "<<FileName<<" position "<<count<<"\n";
				value_checkF="Erreur dans le format du nom de read "+FileName+" position "+str_count;
				exit;
			}
			if(line3.compare("+")==0 and CasavaVersion.compare("1.7")==0){
				control=1;
				cout<<"erreur dans le format du nom de quality seq "<<FileName<<" position "<<count <<" et "<< line3<< "\n";
				value_checkF="Erreur dans le format du nom de quality seq "+FileName+" position "+str_count +" et "+ line3;
				exit;
			}
			if(!strchr("ATGCN",seq[0]))
			{
				control=1;
				cout<<"La sequence "<<FileName<<" position "<<count<<" contient des caracteres non autorisés\n";
				value_checkF="La sequence "+FileName+" position "+str_count+" contient des caracteres non autorisés";
				exit;
			}
			if(line2.size()!=line4.size())
			{
				control=1;
				cout<<"les sequences n'ont pas la même taille dans "<<FileName<<" position "<<count<<"\n";
				value_checkF="les sequences n'ont pas la même taille dans "+FileName+" position "+str_count;
				exit;
			}
			else
			{
				control=0;
				value_checkF="control fichier FASTQ OK";
			}
		}
	}
	fin.clear();
	fin.close();

	cout<<"nb sequences dans " <<read1or2<<"\t"<<(count)<<endl;
	
	if(control==0)
	{
		cout<<"controle de structure fastq du fichier "<<FileName<<" OK\n";
		value_checkC="controle de structure fastq du fichier "+FileName+" OK\n";
	}
	else
	{
		cout<<"verifiez les erreurs dans le fichier "<<FileName<<"\n";
		value_checkC="verifiez les erreurs dans le fichier "+FileName+"\n";	
	}
	
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream checkF;
	checkF<<linenum;	
	timestamp=FormatDate();
	variable=Read+"_Check_FILE";
	value=value_checkF;
	type="TEXT";
	comment="controle contenu de fichier";
	line=checkF.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))
		{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			exit(1);
		}
		else{

		}
	}
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream checkC;
	checkC<<linenum;	
	timestamp=FormatDate();
	variable=Read+"_Check_FILE";
	value=value_checkC;
	type="TEXT";
	comment="controle contenu de fichier";
	line=checkC.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))
		{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			exit(1);
		}
		else{

		}
	}

	char Tcount[20];
	sprintf (Tcount, "%d", count);
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream nbSEQ;
	nbSEQ<<linenum;	
	timestamp=FormatDate();
	variable=Read+"_NBseq";
	value=Tcount;
	type="INTEGER";
	comment="nombre de sequence dans le fichier";
	line=nbSEQ.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	//INSERT INTO DB
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))
		{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			exit(1);
		}
		else{

		}
	}	
	
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream enddb;
	enddb<<linenum;	
	string timestamp_end=FormatDate();
	variable=Read+"_TIME_END";
	value=timestamp_end;
	type="DATETYPE";
	comment="Temps du fin de step";
	line=enddb.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	//INSERT INTO DB
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))
		{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			exit(1);
		}
		else{

		}
	}		
		
	string Fileout=AnalyseREP+"/Verif_"+FileName;	
	fout.open((Fileout).c_str());
	if(fout.good())
	{
		fout<<Fileout<<" "<<control<<endl;
	}
	fout.close();
	return control;	
	mysql_close(connect);	
}

void createStat(string Read, string ID_Processus,map<string,string>SampleData,map<string,string>Connector){
	string line;	
	ifstream fin;
	ofstream fout1;
	ofstream fout2;
	ofstream fout3;
	unsigned long comptread=0;
	Test* membre = new Test;
	string::iterator si;
	map<int,int> nombreNparread;
	int nbNread=0;
	unsigned long nbA=0;
	unsigned long nbT=0;
	unsigned long nbG=0;
	unsigned long nbC=0;
	unsigned long nbN=0;
	unsigned long nbbases=0;
	int index;
	unsigned long nbBases=0;
	string line1;
	string line2;
	string line3;
	string line4;	

	vector<string>ReturnErrorList;	
	string message="";
	string File="";
	string FileName="";
	
	string read1or2="";	
	
	string CasavaVersion="";
	string AnalyseREP="";
	int cyclemax=0;
	string compress="";
	string step="";
	string timestamp="";
	string variable="";
	string value="";
	string type="";
	string comment="";
	string step_number="";
	string query="";
	int linenum=0;
	string repdecomp="";
	string repbrut="";	
	
	string ErrorCode="1";
	string ErrorMessage="";
	MYSQL *connect;
	MYSQL *status;
	MYSQL_ROW row;
	MYSQL_RES *result;
	string readname="";
	size_t found;
	size_t found2;
	string server=Connector["host"];
	string database=Connector["dbname"];
	string user=Connector["user"];
	string pass=Connector["password"];
	string timestamp_deb_stat=FormatDate();

	connect=mysql_init(NULL);
	char* cmd;	
		
		AnalyseREP=SampleData["REPANAL"];
	cyclemax=atoi(SampleData["Cycle"].c_str());
	
	map<string,string>::iterator i ;
	for(i=SampleData.begin();i!=SampleData.end();i++){
		if((*i).first.compare("REPDECOMP")==0){
			repdecomp=(*i).second;
		}
		if((*i).first.compare("REPBRUT")==0){
			repbrut=(*i).second;
		}
	}
	if(repdecomp.compare(repbrut)==0 ){
		if(Read.compare("Read1")==0){
				FileName=SampleData["Anal_Read1"];
								read1or2="R1"; 
								
		}
		if( Read.compare("Read2")==0){
				FileName=SampleData["Anal_Read2"];
								read1or2="R2"; 
								
		}	
	}
	else{
		if(Read.compare("Read1")==0){
			FileName=SampleData["R1_Filename"];
						read1or2="R1"; 
			
		}
		if(Read.compare("Read2")==0){
			FileName=SampleData["R2_Filename"];
						read1or2="R2"; 
						
		}

		CasavaVersion=SampleData["Casava_Version"];
	}
	string Path=repdecomp+"/"+FileName;
	
	int tabscoreread1[cyclemax][65];
	for(int i=0;i<cyclemax;i++){
		for(int j=0;j<65;j++){
			tabscoreread1[i][j]=0;
		}
	}
	string Fileout1 = FileName + "_boxplotvalues.txt";
	string Fileout2 = FileName + "_ATGCN_stats.txt";
	string Fileout3 = FileName + "_Phred_hashtable.txt";

	fin.open(Path.c_str());
	string value_defstat="";
	if(fin.is_open()){
		cout<<"fichier createStat "<<FileName<<" ouvert\n";
		value_defstat="fichier createStat "+FileName+" ouvert";
	}
	else{
		cout<<"probleme à l'ouverture du fichier\n";
		value_defstat="probleme fichier createStat "+Path;
	}
	
	// INSERTION INTO DATABASE //

	connect=mysql_init(NULL);	
	string queryStep="call PC_GetNextStep(@idProcessus:='"+ID_Processus+"');";
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,queryStep.c_str())){
			cout<<"requete envoyee"<<endl;
			exit(1);
		}
		result=mysql_store_result(connect);
              
		if((row=mysql_fetch_row(result))){
			step_number=row[0];
		}
		else{
			cout<<"pas de resultats de requete "<<endl;
		}
	}	
	else{
		ErrorMessage="probleme de connexion a la Database PC_GetNextStep";
		ErrorCode="1";
		cout<<ErrorCode<<endl;
		cout<<ErrorMessage<<endl;
		exit(1);
	}	
	cout<<"step num\t"<< step_number<<endl;

	linenum++;
	stringstream tsart;
	tsart<<linenum;	
	string timestamp_debStat=FormatDate();
	variable=Read+"_TIME_START";
	step="CreateStat";
	timestamp=FormatDate();
	value=timestamp_deb_stat;
	type="DATETYPE";
	comment="Temps du debut du step";
	line=tsart.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	//INSERT INTO DB
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))
		{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			exit(1);
		}
		else{

		}
	}
	//END INSERT INTO DB
	
	linenum++;
	stringstream sif;
	sif<<linenum;	
	timestamp=FormatDate();
	variable=Read+"_FILE";
	step="CreateStat";
	value=value_defstat;
	type="TEXT";
	comment="controle ouverture de fichier read pour createstat";
	line=sif.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	//INSERT INTO DB
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))
		{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			exit(1);
		}
		else{

		}
	}
	//END INSERT INTO DB
	
	while(!fin.eof()){
		getline(fin,line1);
		getline(fin,line2);
		getline(fin,line3);
		getline(fin,line4);
		comptread++;
		index=0;
		nbNread=0;
		for(si=line2.begin();si!=line2.end();si++){
			if(*si=='A'){
				nbA++;
				nbBases++;
			}
			else if(*si=='T'){
					nbT++;
					nbBases++;
			}
			else if(*si=='G'){
					nbG++;
					nbBases++;
			}
			else if(*si=='C'){
					nbC++;
					nbBases++;
			}
			else if(*si=='N'){
					nbN++;
					nbBases++;
					nombreNparread[comptread]++;
			}
			// Calculation of phred score is different in Casava 1.8 and version before
			if(CasavaVersion.compare("1.8")==0)
			{
				int qscore= membre->calcPhredScoreIndexSanger(line4,index);
				tabscoreread1[index][qscore]++;
				index++;
			}
			else
			{
				int qscore= membre->calcPhredScoreIndex(line4,index);
				tabscoreread1[index][qscore]++;
				index++;
			}
		}
	}
	fin.close();
	string value_checkFOUT2="";
	string value_checkFOUT1="";
	string value_checkFOUT3="";
	nbbases=long((comptread-1)*(cyclemax));
	fout2.open((AnalyseREP+'/'+Fileout2).c_str());
	if(fout2.is_open()){
		cout << "Ecriture ATGCN stats "<<read1or2<<"\t"<<Fileout2<<endl;
		value_checkFOUT2="Ecriture ATGCN stats "+read1or2+" "+Fileout2;	
	}
	else{
		cout<<"Erreur d'ouverture "<<Fileout2<<endl;
		value_checkFOUT2="Erreur d'ouverture"+Fileout2;
	}

	// INSERTION INTO DATABASE //
	linenum++;
	stringstream checkF2;
	checkF2<<linenum;	
	timestamp=FormatDate();
	variable=Read+"_Open_FILE";
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
			exit(1);
		}
		else{

		}
	}
			
	// Amount of ATGCN 

	fout2 <<"composition en ATGCN du "<<FileName<<":"<<endl;
	fout2<<"nb total bases calc "<<nbbases<<endl;
	fout2<<"nb total bases  "<<nbBases<<endl;
	fout2<<"nb reads avec N "<<(nombreNparread.size()-1)<<endl;
	fout2 << "#A : "<< nbA<<" "<<(100*((double)nbA/(double)nbbases))<<"%"<< endl;
	fout2 << "#C : "<< nbC<<" "<<(100*((double)nbC/(double)nbbases))<<"%"<< endl;
	fout2 << "#G : "<< nbG<<" "<<(100*((double)nbG/(double)nbbases))<<"%"<< endl;
	fout2 << "#T : "<< nbT<<" "<<(100*((double)nbT/(double)nbbases))<<"%"<< endl;
	fout2 << "#N : "<< nbN<<" "<<(100*((double)nbN/(double)nbbases))<<"%"<< endl;
	
	fout2.close();
	

	fout1.open((AnalyseREP+'/'+Fileout1).c_str());
	if(fout1.is_open()){
		cout << "Ecriture boxplotvalues "<<read1or2<<"\t"<<Fileout1<<endl;
		value_checkFOUT1="Ecriture boxplotvalues "+read1or2+" "+Fileout1;	
	}
	else{
		cout<<"Erreur d'ouverture "<<Fileout1<<endl;
		value_checkFOUT1="Erreur d'ouverture "+Fileout1;
	}
	fout3.open((AnalyseREP+'/'+Fileout3).c_str());
	if(fout3.is_open()){
		cout << "Ecriture phred hashtable "<<read1or2<<"\t"<<Fileout3<<endl;
		value_checkFOUT3="Ecriture phred hashtable "+read1or2+" "+Fileout3;	
	}
	else{
		cout<<"Erreur d'ouverture "<<Fileout3<<endl;
		value_checkFOUT3="Erreur d'ouverture"+Fileout3;
	}
	
	// INSERTION INTO DATABASE //
	linenum++;
	stringstream checkF1;
	checkF1<<linenum;	
	timestamp=FormatDate();
	variable=Read+"_Open_FILE";
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
			exit(1);
		}
		else{

		}
	}

	// INSERTION INTO DATABASE //
	linenum++;
	stringstream checkF3;
	checkF3<<linenum;	
	timestamp=FormatDate();
	variable=Read+"_Open_FILE";
	value=value_checkFOUT3;
	type="TEXT";
	comment="ouverture fichier "+Fileout3;
	line=checkF3.str();
	connect=mysql_init(NULL);
    query="call PC_InsertIntoDBLog('"+ID_Processus+"','"+ step+"','"+step_number+"','"+ timestamp+"','"+ variable+"','"+value+"','"+type+"','"+comment+"','"+line+"');";
	
	if(mysql_real_connect(connect, server.c_str(),user.c_str(), pass.c_str(), database.c_str(), 0, NULL,CLIENT_MULTI_STATEMENTS)){
		if(mysql_query(connect,query.c_str()))
		{
			fprintf(stderr, " %s\n", mysql_error(connect));
			cout<<"probleme de connexion :"+query+"\n";
			exit(1);
		}
		else{

		}
	}		
	
	//creation of boxplot values
	for(int i=0;i<cyclemax;i++)
	{
		int min=0;
		int med=0;
		int quart=0;
		int quart3=0;
		int max=0;
		int sum=0;
		//int count=0;
		long total=0;
		long jphred=0;
		int moyenne=0;
		for(int j=0;j<65;j++){
			int nombre=tabscoreread1[i][j];
			jphred=(long)nombre*j;
			fout3<<i<<" "<<j<<" "<<nombre<<" "<<jphred<<endl;
			sum =nombre+sum;
			total=total+jphred;
			fout3<<i<<" "<<total<<" "<<nombre<<" "<<j<<endl;
			if(nombre!=0){
				if(sum==nombre){min=j;}
				if((int)sum>=(comptread-1)  && quart3==0 && med==0 && quart==0){max=j;quart3=j;med=j;quart=j;						}
				if((int)sum>=(comptread-1)  && med==0 && quart3==0 && quart!=0){max=j;quart3=j;med=j;			}
				if((int)sum>=(comptread-1)  && quart3==0 && med!=0){max=j;quart3=j;			}
				if((int)sum>=(comptread-1)  && quart3!=0 && max==0){max=j;				}
				if((int)sum<(comptread-1) && (int)sum>=(comptread-1)*3/4 && quart3==0 && med==0 && quart==0){quart3=j;med=j;quart=j;				}
				if((int)sum<(comptread-1) && (int)sum>=(comptread-1)*3/4 && med!=0 && quart3==0){quart3=j;				}
				if((int)sum<(comptread-1) && (int)sum>=(comptread-1)*3/4 && med==0 && quart3==0 && quart!=0){quart3=j;med=j;				}
				if((int)sum<(comptread-1)*3/4 && (int)sum>=(comptread-1)/2 && quart==0 && med==0){med=j;quart=j;}
				if((int)sum<(comptread-1)*3/4 && (int)sum>=(comptread-1)/2 && quart!=0 && med==0){med=j;}
				if((int)sum<(comptread-1)/2 && (int)sum>=(comptread-1)/4 && quart==0){quart=j;}
			}
		}
		moyenne=(int)(total/(comptread-1));
		fout1<<i+1<<" "<<min<<" "<<quart<<" "<<med<<" "<<quart3<<" "<<max<<" "<<moyenne<<endl;
		fout3<<"\n";
	}

	fout3.close();
	fout1.close();

	// INSERTION INTO DATABASE //
	linenum++;
	stringstream endT;
	endT<<linenum;	
	string timestamp_end=FormatDate();
	cout<<"FIN DU STEP\t"<< timestamp_end<<endl;
	variable=Read+"_TIME_END";
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
			exit(1);
		}
		else{

		}
	}	
	mysql_close(connect);
}
