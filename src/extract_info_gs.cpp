#include<iostream>
#include<fstream>
#include<sstream>
#include<stdlib.h>

using namespace std;

int main(int argc,char *argv[])
{
//	int linenum=atoi(argv[2]);
	ifstream readfile(argv[1],ios::in);
	if(!readfile)
	cerr<<"read the file failed"<<endl;
	string line;
	int k;
	k=0;
	while(!readfile.eof())
	{
		getline(readfile,line);
		if(line.size()==0)
		continue;
		stringstream ss;
		ss<<line;
		int num;
		ss>>num;
		if(k>num){
			break;
		}
		cout<<num<<endl;
		k=num;
	}
	return 0;
}
