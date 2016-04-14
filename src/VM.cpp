// Implementation of Vertex Mover Algorithm
// WARNING: The network is assumed to be undirected
// edge file has to be not symmetrized
// Simultaneous presence of both links (a b weigth and b a weight) causes 
// wrong program output

#include<iostream>
#include<set>
#include<vector>
#include<fstream>
#include<string>
#include<sstream>
#include <stdlib.h>
#include <algorithm>

using namespace std;

struct link2{
  unsigned int partner;
  double value;
  
  link2(){
  };

  link2(const double & Q, const unsigned int & part){
    this->value=Q;
    this->partner = part;
  };
};

struct partner_sorted{
  bool operator()(const link2 & a, const link2 & b){
    return (a.partner > b.partner);
  };
};

struct value_sorted{
 bool operator()(const link2 & a, const link2 & b){
    return ((a.value<b.value)||((a.value==b.value)&&(a.partner < b.partner)));
 };
};


class NO{
public:
  int read_in(string nodefilename,string linkfilename);
  void print_link_matrix();
  void print_order();
  void perform_NO(string output_file);
  
private:
  double Q;
  double total_link_weight;

  vector<double> indegree;
  vector<double> outdegree;
  vector<double> degree;
  
  vector<double> community_degree;
  vector<double> community_indegree;
  vector<double> community_outdegree;

  vector<unsigned int> colors;
  
  bool algorithm_round();
  bool algorithm_init();

  vector<set<link2,partner_sorted> > link_matrix;
  set<link2,value_sorted> order;


};

int NO::read_in(string nodefilename,string linkfilename){
  fstream hFile(linkfilename.c_str());
  
  if (!hFile){
    cerr << "File " << linkfilename << " could not be opened. Will exit!"<< endl;
    exit(0);
  }

 fstream hNode(nodefilename.c_str());
  
  if (!hNode){
    cerr << "File " << nodefilename << " could not be opened . Will exit!" << endl;
    exit(0);
  }

  unsigned int node = 0;
  unsigned int color=0;
  unsigned int dmp=0;
  
  // Empty color = 0

  while (!hNode.eof()){
    hNode >> node >> color >> dmp;
    if (node >= colors.size()){
      colors.resize(node+1,0);
    };
    
    colors[node]=color+1;
  };


  // Initiate degree sorted list

  const unsigned int color_size=colors.size();
  indegree.resize(color_size,0);
  outdegree.resize(color_size,0);
  degree.resize(color_size,0);



  int start=-1;
  int ende=-1;
  double value=-1;

  {
     set<link2,partner_sorted> empty_set;
     link_matrix.resize(color_size,empty_set);
  }

  while (!hFile.eof()){
    hFile >> start >> ende >> value;
    if (value != 0.0){
      int max=start<ende?ende:start;
      int min=start<ende?start:ende;
      // only valid for symmetric implementation
      total_link_weight+=2*value;

      if (link_matrix.size()<=max){
	cerr << "Links to not existing node. Will exit" << endl;
	exit(0);
      }

      link2 element_to_insert(value,min);
      link_matrix[max].insert(element_to_insert);
      element_to_insert.partner=max;
      link_matrix[min].insert(element_to_insert);
      outdegree[start]+=value;
      indegree[ende]+=value;
      degree[start]+=2*value;
      degree[ende]+=2*value;
      
    }
  }

  // compensating the double read of the last line 
  total_link_weight-=2*value;
  outdegree[start]-=value;
  indegree[ende]-=value;
  degree[start]-=2*value;
  degree[ende]-=2*value;


  for (unsigned int inode=0;inode < color_size;inode++){
    link2 to_insert(degree[inode],inode);
    order.insert(to_insert);
  }

  //test
  /*  for (unsigned int iruni=0;iruni< degree.size();iruni++){
    cout << iruni << "\t" << degree[iruni] << endl;
    }*/

  return 0 ;

  

}

void NO::print_order(){
  cout << "The order of the elements I would consider is" << endl;
  unsigned int counter=0;
  for (set<link2,value_sorted>::iterator iter = order.begin();iter != order.end();iter++){
    cout << counter++ << "\t Node " << iter->partner << "\t degree : " << iter->value << endl;
  }

}

void NO::print_link_matrix(){
  const unsigned int sizer=link_matrix.size();
  for (unsigned int inode=0;inode < sizer; inode ++){
    cout << endl << inode << ":"; 

    for (set<link2,partner_sorted>::iterator iter=link_matrix[inode].begin();iter != link_matrix[inode].end();iter++){

      cout << iter->partner << " (" << iter->value << " )\t" ;

    };
  };

  cout << endl;

};

bool NO::algorithm_init(){
  vector <double> internal_links;
  const unsigned int sizer=colors.size();
  for (unsigned int inode=0;inode < sizer;inode ++){
    const unsigned int color=colors[inode];
    if (internal_links.size()<=color){
      internal_links.resize(color+1,0);
      community_degree.resize(color+1,0);
      community_indegree.resize(color+1,0);
      community_outdegree.resize(color+1,0);
    }
    
    for (set<link2,partner_sorted>::iterator iter=link_matrix[inode].begin();iter!=link_matrix[inode].end();iter++){
	const unsigned int partner=iter->partner;
	const double value=iter->value;
	const unsigned int partner_color=colors[partner];
      if (partner_color==color){
	internal_links[color]+=value;
      }
      
      if (internal_links.size()<=partner_color){
	internal_links.resize(partner_color+1,0);
	community_degree.resize(partner_color+1,0);
	community_indegree.resize(partner_color+1,0);
	community_outdegree.resize(partner_color+1,0);
      }
      


      community_indegree[partner_color]+=value;
      community_outdegree[color]+=value;
      community_degree[partner_color]+=value;
      community_degree[color]+=value;
    };
  }
  
  Q=0.0;
  const unsigned col_sizer=community_degree.size();
  for (unsigned int icolor=0;icolor<col_sizer;icolor++){
    // Only for symmetric variant 
    Q+=internal_links[icolor]/total_link_weight - community_degree[icolor]/total_link_weight * community_degree[icolor]/total_link_weight/4.0;
  }
};

bool NO::algorithm_round(){
  bool node_changed=false;
  double Qchange=0.0;

  for (set<link2,value_sorted>::iterator iter=order.begin();iter!=order.end();iter++){
    const unsigned int inode=iter->partner;

    //Only valid for symmetric from here on
    vector<unsigned int> partner_colors;
    
    vector<double> links_from_node;

    for (set<link2,partner_sorted>::iterator iterPartner=link_matrix[inode].begin();iterPartner!=link_matrix[inode].end();iterPartner++){
      const unsigned int partner_color=colors[iterPartner->partner];
      const double value=iterPartner->value;

      if (partner_color>=links_from_node.size()){
	links_from_node.resize(partner_color+1,0);
      }

      if (links_from_node[partner_color]==0){
	partner_colors.push_back(partner_color);
      }

      links_from_node[partner_color]+=2*value;

    }


    //to obtain consistency
    sort(partner_colors.begin(),partner_colors.end());


    double Qmax=-1e19;

    unsigned int opt_color=0;

    unsigned int max_partner=partner_colors.size();
   
    const unsigned int home_color=colors[inode];
 
    const double links_from_home=links_from_node.size()<=home_color?0:links_from_node[home_color];

    const double degree_home_community=community_degree[home_color];
   
    for (unsigned int ipartner=0; ipartner< max_partner;ipartner++){
      const unsigned int partner_color=partner_colors[ipartner];

      double tempQ=partner_color!=home_color?((links_from_node[partner_color] - links_from_home)/total_link_weight 
	- degree[inode]/2.0/total_link_weight*(community_degree[partner_color]-degree_home_community+ degree[inode])/total_link_weight):0;
      
      if (Qmax<tempQ){
	opt_color=partner_color;
	Qmax=tempQ;
      }
      
    }

    if ((opt_color!=0)&&(Qmax>0.0)){
      Qchange+=Qmax;
       cout << "for node " << inode << " optimal value " << Qmax
	    << " from color " << colors[inode] <<  " with color " << opt_color << endl;

      node_changed=true;
      community_degree[colors[inode]]-=degree[inode];
      colors[inode]=opt_color;
      community_degree[opt_color]+=degree[inode];
    }

    const unsigned int sizer=partner_colors.size();
    //Put it back in original state
    for (unsigned int inode=0;inode<sizer;inode++){
      links_from_node[partner_colors[inode]]=0;

    }
    
    partner_colors.clear();
  }
  Q+=Qchange;

  return node_changed;
};

void NO::perform_NO(string output_file){
  algorithm_init();
  cout << "Initial Q  " << Q << endl;

  unsigned int counter=0;
  while (algorithm_round()){cout << "Runde " << counter++ << endl;};
  cout << "Final solution after " << counter+1 << " rounds with Q= " << Q << endl;
  const unsigned int sizer=colors.size();
      stringstream filename;
      //filename << output_file << "-boosted-to-" << Q;
			filename<< output_file;
   ofstream file(filename.str().c_str());
	
	file << Q <<endl;
  for (unsigned int inode=0;inode < sizer;inode++){
    file << inode << "\t" << colors[inode] << endl;
    }
};

int main(int argc,char* argv[]){
  string nodefilename(argv[1]);
  string linkfilename(argv[2]);
  NO no;
  cout << "Will read in " << nodefilename << " (nodes)  and " 
       << linkfilename << " (links) with status " 
       << no.read_in(nodefilename,linkfilename) <<  endl;
	
	cout<<nodefilename<<endl;
  string output_file="temp_MSG_result_of_";
	string temp;
	temp.assign(nodefilename,8,nodefilename.length()) ;
	output_file += temp;
	//cout<<output_file<<endl;
  no.perform_NO(output_file);
}
