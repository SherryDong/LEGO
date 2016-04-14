// FINAL VERSION: USED TO COMPILE STATIC VARIANT FROM January 26,2007
// // ATTENTION: Results differe if compiled with 64 bit or 32 bit support
// // The reference variant was always compiled with 32 bit support only
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <map>
#include <vector>
#include <string>
#include <stdlib.h>
#include <algorithm>

using namespace std;

#define DEBUG
#define VM_dir "/home/sherry/project/network/NETGO/LEGO_code_v6.0_RNA/bin//"

const long MAX_NUMBER_OF_NODE = 30000;

int N = 0;

map<string, int> name2node;
map<int, string> node2name;

vector< vector<string> > answer;

float connection[MAX_NUMBER_OF_NODE][MAX_NUMBER_OF_NODE] = {0};

string input_file;

void dealWithMatrix(string matrix_file)
{
	ifstream in(matrix_file.c_str());
	if(!in)
	{
		cout<<"matrix file does not exist."<<endl;
		exit(-1);
	}

	int count = 0;
	string start, end;
	float value;
	int index = 0;
	int index2 = 0;

	ofstream out_matrix(("matrix_"+matrix_file).c_str());
	ofstream out_vm_matrix(("vm_matrix_"+matrix_file).c_str());

	while(in >> start >> end >> value)
	{
		if(start == end)
		{
			++ count;
			continue;
		}
		map<string, int>::iterator temp = name2node.find(start);
		if(temp == name2node.end())
		{
			name2node.insert(make_pair(start, ++index));
			node2name.insert(make_pair(++index2, start));
		}
		temp = name2node.find( end);
		if(temp == name2node.end())
		{
			name2node.insert(make_pair(end, ++index));
			node2name.insert(make_pair(++index2, end));
		}
		//out_matrix<<name2node[start]<<'\t'<<name2node[end]<<'\t'<<value<<endl;
		//out_matrix<<name2node[end]<<'\t'<<name2node[start]<<'\t'<<value<<endl;
		//out_vm_matrix<<name2node[start]<<'\t'<<name2node[end]<<'\t'<<value<<endl;

		connection[name2node[start]][name2node[end]] = value;
		connection[name2node[end]][name2node[start]] = value;
	}
	for(map<string, int>::iterator it_i=name2node.begin(); it_i!=name2node.end(); ++it_i)
		if(it_i->second > N)
			N = it_i->second;
	for(int i=1; i<N; ++i)
		for(int j=i+1; j<=N; ++j)
			if(connection[i][j] > 0)
			{
				out_matrix<<i<<'\t'<<j<<'\t'<<connection[i][j]<<endl;
				out_matrix<<j<<'\t'<<i<<'\t'<<connection[i][j]<<endl;
				out_vm_matrix<<i<<'\t'<<j<<'\t'<<connection[i][j]<<endl;

			}

	ofstream out_node(("node_"+matrix_file).c_str());
	for(map<int, string>::iterator it_i=node2name.begin(); it_i!=node2name.end(); ++it_i)
		out_node<<it_i->first<<'\t'<<it_i->second<<endl;
}

struct data
{
	double deltaQ;
	unsigned int partner;

	data(){}

	data(const double & delta, const unsigned int & part)
	{
		this->deltaQ=delta;
		this->partner = part;
	}
};

struct data_all
{
	double deltaQ;
	unsigned int start;
	unsigned int ende;

	data_all(const double & delta, const unsigned int & start, const unsigned int & ende)
	{
		this->deltaQ=delta;
		this->start=start;
		this->ende=ende;
	}
};


struct data_min
{
	double deltaQ;
	unsigned int partner;
	bool min;

	data_min(){}

	data_min(const double & delta, const unsigned int & part,const bool & min)
	{
		this->deltaQ=delta;
		this->partner = part;
		this->min= min;
	}
};

struct compareme
{

	bool operator()(const  data & a, const data & b)
	{
		return ((a.deltaQ > b.deltaQ)||((a.deltaQ==b.deltaQ)&&(a.partner<b.partner)));
	}
};

struct compareme_index
{
	 bool operator()(const data_min & a, const data_min &b)
	 {
		 return ((a.partner < b.partner)||((a.partner==b.partner)&&(a.min<b.min)));
	 }
};

struct compare_pairs
{
	bool operator()(const data_all & a, const data_all & b)
	{
		return ((a.deltaQ>b.deltaQ)||((a.deltaQ==b.deltaQ)&&(a.start<b.start))||((a.deltaQ==b.deltaQ)&&(a.start==b.start)&&(a.ende< b.ende)));
	}
};


struct compare_pairs_no_deltaQ
{
	bool operator()(const data_all & a, const data_all & b)
	{
		return ((a.start<b.start)||((a.start==b.start)&&(a.ende< b.ende)));
	}
};



class community
{
public:
	community(){init(1);}

	void init(unsigned int nodes);
	int print_colors(string outputfilename);
	 void merge(unsigned int comm1,unsigned int comm2);
	 ~community(){communities.clear();}

private:
	vector<vector<unsigned int> > communities;

};

void community::init(unsigned int nodes)
{
	communities.clear();
	for (unsigned int inode=0;inode<nodes;inode++)
	{
		 vector<unsigned int> initial = vector<unsigned int>(1,inode);
		 communities.push_back(initial);
	}
}

void community::merge(unsigned int comm1, unsigned int comm2)
{
	const unsigned int sizer=communities[comm2].size();
	for (unsigned int inode=0;inode<sizer;inode++)
		communities[comm1].push_back(communities[comm2][inode]);
	communities[comm2].clear();
}

int community::print_colors(string outputfilename)
{
	unsigned int counter=0;
	const unsigned int sizer=communities.size();
	vector<unsigned int> colors(sizer,0);
	for (unsigned int icomm=0;icomm<sizer;icomm++)
	{
		const unsigned int sizer2=communities[icomm].size();
		for (unsigned int inode=0;inode<sizer2;inode++)
			colors[communities[icomm][inode]]=counter;
		if (sizer2>0)
			counter++;
	}

   cout << " Will output into " << outputfilename << endl;

   ofstream file(outputfilename.c_str());
	if (!file)
	{
		cerr << " Something went wrong with the output. Will exit " << endl;
		exit(0);
	}

	const unsigned int sizer3=colors.size();
	for (unsigned int inode=0;inode<sizer3;inode++)
		file << inode << "\t" << colors[inode] << "\t1" << endl;
	cout << "Number of colors " << counter << endl;
}








class algorithm
{
public:
	double do_CMS();
	double do_best_merges(const string &  filename, const unsigned int & levels);
	int read_in(string filename);



private:
	void merge_two_colors(const unsigned int & color_1, const unsigned int & color_2);
	vector<set<data,compareme> > deltaQ;
	vector<double> indegree;
	vector <double> outdegree;
	vector <double> degree;
	double total_link_weight;
	community com;
};


void algorithm::merge_two_colors(const unsigned int & start, const unsigned int & ende)
{
	const unsigned int min=start<ende?start:ende;

	const unsigned int max=start<ende?ende:start;

	// Merging two communities in community list
	com.merge(min,max);

	set<data_min,compareme_index> merging_structure;
	for (set<data,compareme>::iterator iter=deltaQ[min].begin();iter!=deltaQ[min].end();iter++)
	{
		if (iter->partner!=max)
		{
			data_min dmmy(iter->deltaQ,iter->partner,true);
			merging_structure.insert(dmmy);
		}
	}

	for (set<data,compareme>::iterator iter=deltaQ[max].begin();iter!=deltaQ[max].end();iter++)
	{
		if (iter->partner != min)
		{
			data_min dmmy(iter->deltaQ,iter->partner,false);
			merging_structure.insert(dmmy);

		}
	}

	set<data,compareme> merged_structure;
	{
		set<data_min,compareme_index>::iterator iter=merging_structure.begin();
		set<data_min,compareme_index>::iterator next;
		while(iter!=merging_structure.end())
		{
			data new_element(iter->deltaQ,iter->partner);
			next=iter;
			next++;

			data element_to_delete(iter->deltaQ,iter->min==true?min:max);
			deltaQ[iter->partner].erase(element_to_delete);

			if (iter->partner==next->partner)
			{
				new_element.deltaQ=iter->deltaQ+next->deltaQ;
				//min is the first element
				element_to_delete.partner=max;
				deltaQ[iter->partner].erase(element_to_delete);

				//Removing second element
				element_to_delete.deltaQ=next->deltaQ;
				element_to_delete.partner=min;
				deltaQ[iter->partner].erase(element_to_delete);

				merging_structure.erase(next);
			}
			else
			{
				if (iter->min==true)
					new_element.deltaQ-=(indegree[max]*outdegree[iter->partner]+indegree[iter->partner]*outdegree[max])/total_link_weight/total_link_weight;
				else
					new_element.deltaQ-=(indegree[min]*outdegree[iter->partner]+indegree[iter->partner]*outdegree[min])/total_link_weight/total_link_weight;
			}


			if (iter->min==false)
			{
				data element_to_delete(iter->deltaQ,min);
				deltaQ[max].erase(element_to_delete);
				element_to_delete.deltaQ=new_element.deltaQ;
				deltaQ[max].insert(element_to_delete);
			}


			merged_structure.insert(new_element);

			// Remove the remanent colors


			element_to_delete.partner=min;
			element_to_delete.deltaQ=new_element.deltaQ;

			deltaQ[iter->partner].insert(element_to_delete);

			++iter;
		}

	}

	//Delete bigger entry and add smaller entry
	deltaQ[max].clear();
	deltaQ[min]=merged_structure;

	indegree[min]+=indegree[max];
	outdegree[min]+=outdegree[max];
	indegree[max]=0;
	outdegree[max]=0;

	degree[min]+=degree[max];
	degree[max]=0;
}

int algorithm::read_in(string filename)
{
	fstream file(filename.c_str());

	int start=-1;
	int ende=-1;
	unsigned int global_max=0;
	double value=-1.0;
	total_link_weight=0.0;

	while (! file.eof())
	{
		file >> start >> ende >> value;
		int max=start>ende?start:ende;

		if (indegree.size() <= max)
		{
			indegree.resize(max+1,0.0);

			outdegree.resize(max+1,0.0);

			degree.resize(max+1,0.0);
		}

		if (!file.eof())
		{
			indegree[ende]+=value;
			outdegree[start]+=value;
			degree[start]+=value;
			degree[ende]+=value;
			if (max>global_max)
				global_max=max;

			total_link_weight+=value;
		}
	}

  //  cout << "Total Link weight " << total_link_weight << endl;

	file.close();

	file.open(filename.c_str());

	deltaQ.clear();

	set<data,compareme> dmmy;

	deltaQ.resize(global_max+1,dmmy);

	cout << "Warning I assume that all link definition only appear once. Otherwise I am lost " << endl;
	data new_element;

	while (! file.eof())
	{
		file >> start >> ende >> value;

		// Search whether something is present
		unsigned to_consider=deltaQ[start].size()>deltaQ[ende].size()?ende:start;
		unsigned to_search=to_consider==ende?start:ende;

		if (!file.eof())
		{
			double found_costfunction=0.0;
			bool found = false ;
			set<data,compareme>::iterator iter=deltaQ[to_consider].begin();
			while (iter!=deltaQ[to_consider].end())
			{
				if (iter->partner==to_search)
				{
					found_costfunction=iter->deltaQ;
					iter++;
					deltaQ[to_consider].erase(data(found_costfunction,to_search));
					data element_to_delete(found_costfunction,to_consider);
					deltaQ[to_search].erase(element_to_delete);
					found = true;
				}
				else
					iter++;
			}


			if (found)
				new_element.deltaQ=found_costfunction+value/total_link_weight;
			else
				new_element.deltaQ=found_costfunction+value/total_link_weight - (indegree[start]*outdegree[ende]+indegree[ende]*outdegree[start])/total_link_weight/total_link_weight;

			new_element.partner=ende;

			deltaQ[start].insert(new_element);

			new_element.partner=start;

			deltaQ[ende].insert(new_element);
		}
	}

  com.init(global_max+1);
  return 0;
}







double algorithm::do_CMS()
{
	const unsigned int number_nodes=degree.size();

	double Q=0.0;

	for (unsigned inode =0; inode<number_nodes;inode++)
		Q-=indegree[inode]/total_link_weight*outdegree[inode]/total_link_weight;

	cout << "Initial Costfunction " << Q << endl;

	unsigned int counter=0;
	double max_deltaQ=-1;

	do
	{
		max_deltaQ=-1;

		unsigned int first=-1;
		unsigned int partner=-1;

		const unsigned int sizer=deltaQ.size();
		for (unsigned int icom=0;icom < sizer;icom++)
		{
			if (deltaQ[icom].size()>0)
			{
				set<data,compareme>::iterator iter=deltaQ[icom].begin();
				if (max_deltaQ < iter->deltaQ)
				{
					first=icom;
					max_deltaQ=iter->deltaQ;
					partner=iter->partner;
				}
			}
		}

		if (max_deltaQ>=0)
		{
			cout << "Merging " << first << " with " << partner << endl;

			merge_two_colors(first,partner);
			Q+=max_deltaQ;
		}

		cout << "Round  " << counter++ << " Q = " << Q << endl;

		}while (max_deltaQ >=0);

  return Q;

};

double algorithm::do_best_merges(const string & filename, const unsigned int & levels)
{
	const unsigned int number_nodes=degree.size();

	const double square_minimal_costfunction_difference=1.0/total_link_weight/total_link_weight/total_link_weight/total_link_weight;

	double maxQfound=-1e19;

	double Q=0.0;

	for (unsigned inode =0; inode<number_nodes;inode++)
		Q-=indegree[inode]/total_link_weight*outdegree[inode]/total_link_weight;


	cout << "Initial Costfunction " << Q << endl;

	unsigned int counter=0;
	double max_deltaQ=0.0;

	set<data_all,compare_pairs> data_container;


	const unsigned int sizer=deltaQ.size();
	for (unsigned int icom=0;icom < sizer;icom++)
	{
		if (deltaQ[icom].size()>0)
		{
			for (set<data,compareme>::iterator iter=deltaQ[icom].begin();iter!=deltaQ[icom].end();iter++)
			{
				if (icom < iter->partner)
				{
					data_all element_to_insert(iter->deltaQ,icom, iter->partner);
					data_container.insert(element_to_insert);
				}
			}
		}
	}


	do
	{
		max_deltaQ=0.0;

		unsigned int first=-1;
		unsigned int partner=-1;

		unsigned int level_counter=1;

		vector<bool> touched(indegree.size(),false);

		double last_costfunction=data_container.begin()->deltaQ;

		set<data_all,compare_pairs> elements_to_delete;
		set<data_all,compare_pairs_no_deltaQ> elements_to_add;

		vector <unsigned int> touched_communities;



		for (set <data_all, compare_pairs>::iterator iter = data_container.begin(); ((iter!=data_container.end())&&(level_counter<=levels)); iter++)
		{
			const unsigned int start=iter->start;
			const unsigned int ende=iter->ende;
			const double temp_deltaQ=iter->deltaQ;

			// Updating the level counter
			if ((last_costfunction-temp_deltaQ)*(last_costfunction-temp_deltaQ) >  square_minimal_costfunction_difference)
			{
				level_counter++;
				last_costfunction=temp_deltaQ;
			}


			if ((level_counter<=levels)&&(!touched[start])&&(!touched[ende])&&(last_costfunction>0))
			{
				set<data,compareme>::const_iterator delta=deltaQ[start].begin();

				while (delta != deltaQ[start].end())
				{
					unsigned int min=start<delta->partner?start:delta->partner;
					unsigned int max=start==min?delta->partner:start;

					data_all element_to_delete(delta->deltaQ,min,max);
					elements_to_delete.insert(element_to_delete);
					delta++;
				}

				for (set<data,compareme>::iterator delta=deltaQ[ende].begin();delta!=deltaQ[ende].end();delta++)
				{
					unsigned int min=delta->partner<ende?delta->partner:ende;
					unsigned int max=ende==min?delta->partner:ende;
					data_all element_to_delete(delta->deltaQ,min,max);
					elements_to_delete.insert(element_to_delete);
				}

				cout << "Round  " << counter << " merging " << start << "\t" << ende << endl;
				merge_two_colors(start,ende);

				touched_communities.push_back(start);
				touched_communities.push_back(ende);

				// add the changed elements again

				Q+=temp_deltaQ;
				max_deltaQ++;
				touched[start]=true;
				touched[ende]=true;
			}
		}


    // Applying changes of deltaQ on datacontainer
		{
			set<data_all,compare_pairs>::iterator ielement=elements_to_delete.begin();
			while (ielement!=elements_to_delete.end())
			{
				data_container.erase(data_all(ielement->deltaQ,ielement->start,ielement->ende));
				ielement++;
			}
		}

    //    cout << "data deleted " << endl;



    //Adding the new data
		{
			const unsigned int sizer = touched_communities.size();
			for (unsigned int icomm=0;icomm < sizer;icomm++)
			{
				const unsigned int start=touched_communities[icomm];
				for (set<data,compareme>::iterator some=deltaQ[start].begin();some != deltaQ[start].end();some++)
				{
					unsigned int min=some->partner;
					unsigned int max=min<start?start:min;
					min=max==min?start:min;

					data_container.insert(data_all(some->deltaQ,min,max));
				}
			}
		}

		cout << "Round " << counter++ << " Q = " << Q << endl;

		if (Q > maxQfound)
			maxQfound=Q;

  }while (max_deltaQ > 0);

  cout << "Best overall solution found " << maxQfound << endl;

  stringstream filestream;

	filestream << "partion_"<<input_file;
  //filestream << "partition-MSG-" << filename << "-l-eq-"<< levels << "-with-Q-"<< maxQfound ;

  cout << "OUTPUTTING PARTITION TO " << filestream.str() << endl;

  com.print_colors(filestream.str());

  return Q;

}

void changeAnswer()
{
	string outfile = "MSG_result_of_" + input_file ;
	ofstream out(outfile.c_str());

	string temp_answer_file = "temp_" + outfile;
	ifstream in(temp_answer_file.c_str());
	if(!in)
	{
		cout<<"temp answer input error"<<endl;
		exit(-1);
	}
	double M = 0;
	in >> M;
	
	string s;
	getline(in,s);
	getline(in,s);

	vector<int> node ;
	vector<int> cluster;


	int temp1, temp2;
	while(in >> temp1 >> temp2)
	{
		//cout<<temp1<<'\t'<<temp2<<endl;
		-- temp2;
		node.push_back(temp1);
		cluster.push_back(temp2);
	}
	vector<int>::iterator it = max_element(cluster.begin(), cluster.end());
	answer.resize(*it);
	
	for(int i=0; i<N; ++i)
		answer[cluster[i]-1].push_back(node2name[node[i]]);

	for(int i=0; i<answer.size(); ++i)
	{
		for(int j=0; j<answer[i].size(); ++j)
			out<<answer[i][j]<<'\t';
		out<<endl;
	}

	out<<M<<'\t'<<endl;
	remove (temp_answer_file.c_str());
//	system(("rm -rf " + temp_answer_file).c_str());
}

int main(int argc,char* argv[])
{
	input_file = argv[1];
	dealWithMatrix(input_file);

	algorithm algo;
	string matrix_file = "matrix_" + input_file;
	string node_file = "node_" + input_file;
	string vm_matrix_file = "vm_matrix_"+input_file; 
	string partion_file = "partion_"+input_file; 
	string vm_out = VM_dir; 
		vm_out += "VM.out"; 

	cout << "Will read in " << input_file << endl;
	algo.read_in(matrix_file);
	unsigned int levels = (int) (0.35 * N);
	if(levels == 0)
		levels = 1;
	cout << "Will perform best_merges with l= " << levels << endl;
	algo.do_best_merges(matrix_file,levels);

	system((vm_out + " "+partion_file+ " "+vm_matrix_file).c_str());
	cout << "finish VM.out running"<< endl; 
	remove(node_file.c_str()); 
	remove(matrix_file.c_str()); 
	remove(vm_matrix_file.c_str()); 
	remove(partion_file.c_str()); 
	//system(("rm -rf matrix_"+input_file+" node_"+input_file+" partion_"+input_file+" vm_matrix_"+input_file).c_str());

	changeAnswer();
}
