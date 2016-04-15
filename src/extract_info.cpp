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
	while(!readfile.eof())
	{
		getline(readfile,line);
		if(line.size()==0)
		continue;
		stringstream ss;
		ss<<line;
		int num;
		ss>>num;
		double a[3];
		for(int i=0;i<3;i++)
		ss>>a[i];
		cout<<a[2]<<endl;
	}
	return 0;
}
