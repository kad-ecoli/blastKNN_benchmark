const char* docstring="\n"
"calculate_ic uniprot_exp.F.is_a naive.F is_a.csv name.csv\n"
"    calculate the naive probability and information content of GO terms\n"
"\n"
"Input:\n"
"    uniprot_exp.F.is_a - GO annotation in the following format\n"
"                         [1]  accession\n"
"                         [2]  comma separated list of GO terms, with parent terms\n"
"    is_a.csv           - GO_ID Aspect is_a_direct is_a_indirect\n"
"    name.csv           - GO_ID Aspect name\n"
"\n"
"Output:\n"
"    naive.F            - naive probability and information content\n"
"                         [1] GO_ID\n"
"                         [2] naive probability\n"
"                         [3] information content\n"
"                         [4] information content condition on parents\n"
;

#include <iostream>
#include "StringTools.hpp"

using namespace std;

void parse_is_a_file(const string &is_a_filename, 
    map<string,vector<string> >&is_a_dict)
{
    bool fromStdin=(is_a_filename=="-");
    ifstream fin;
    vector<string> line_vec;
    string line;
    string GOterm,is_a_direct;
    int i;
    if (!fromStdin) fin.open(is_a_filename.c_str());
    while((fromStdin)?cin.good():fin.good())
    {
        if (fromStdin) getline(cin,line);
        else           getline(fin,line);
        if (line.size()==0) continue;
        Split(line, line_vec, '\t');
        GOterm       =line_vec[0];
        is_a_direct  =line_vec[2];
        for (i=0;i<line_vec.size();i++) line_vec[i].clear(); line_vec.clear();
        Split(is_a_direct,line_vec,',');
        is_a_dict[GOterm]=line_vec;
        for (i=0;i<line_vec.size();i++) line_vec[i].clear(); line_vec.clear();
        GOterm.clear();
        is_a_direct.clear();
    }
    if (!fromStdin) fin.close();
    vector<string> ().swap(line_vec);
    string ().swap(line);
    string ().swap(GOterm);
    string ().swap(is_a_direct);
}

void parse_name_file(const string &name_filename, 
    map<string,string>&name_dict)
{
    bool fromStdin=(name_filename=="-");
    ifstream fin;
    vector<string> line_vec;
    string line;
    string GOterm,name;
    int i;
    if (!fromStdin) fin.open(name_filename.c_str());
    while((fromStdin)?cin.good():fin.good())
    {
        if (fromStdin) getline(cin,line);
        else           getline(fin,line);
        if (line.size()==0) continue;
        Split(line, line_vec, '\t');
        GOterm=line_vec[0];
        name  =line_vec[2];
        for (i=0;i<line_vec.size();i++) line_vec[i].clear(); line_vec.clear();
        name_dict[GOterm]=name;
        GOterm.clear();
        name.clear();
    }
    if (!fromStdin) fin.close();
    vector<string> ().swap(line_vec);
    string ().swap(line);
    string ().swap(GOterm);
    string ().swap(name);
}


void calculate_ic(const string &inputfilename,
    const string &outputfilename,
    map<string,vector<string> > &is_a_dict,
    map<string,string> &name_dict)
{
    map<string,vector<string> > allterm_dict;
    bool fromStdin=(inputfilename=="-");
    ifstream fin;
    vector<string> line_vec;
    vector<string> accession_list;
    string line;
    string accession,GOterm,GOterms;
    size_t i;
    if (!fromStdin) fin.open(inputfilename.c_str());
    while((fromStdin)?cin.good():fin.good())
    {
        if (fromStdin) getline(cin,line);
        else           getline(fin,line);
        if (line.size()==0) continue;
        Split(line, line_vec, '\t');
        accession=line_vec[0];
        GOterms  =line_vec[1];
        for (i=0;i<line_vec.size();i++) line_vec[i].clear(); line_vec.clear();
        Split(GOterms, line_vec, ',');
        allterm_dict[accession]=line_vec;
        for (i=0;i<line_vec.size();i++) line_vec[i].clear(); line_vec.clear();
        accession_list.push_back(accession);
    }
    if (!fromStdin) fin.close();
    vector<string> ().swap(line_vec);

    vector<string> GOterm_list;
    map<string,vector<string> > GOterm_dict;
    vector<string> tmp_vec;
    size_t a;
    for (a=0;a<accession_list.size();a++)
    {
        accession=accession_list[a];
        for (i=0;i<allterm_dict[accession].size();i++)
        {
            GOterm=allterm_dict[accession][i];
            if (GOterm_dict.count(GOterm)==0)
            {
                GOterm_list.push_back(GOterm);
                GOterm_dict[GOterm]=tmp_vec;
            }
            GOterm_dict[GOterm].push_back(accession);
        }
    }
    vector<string> ().swap(accession_list);

    vector<pair<double,string> > prob_vec;
    map<string,double> ic_dict;
    map<string,double> ic_condition_dict;
    map<string,size_t> accession_count;
    size_t g;
    double prob;
    size_t num_with_parent=0;
    string parent;
    for (g=0;g<GOterm_list.size();g++)
    {
        GOterm=GOterm_list[g];
        prob=(1.0+GOterm_dict[GOterm].size())/(1.0+allterm_dict.size());
        prob_vec.push_back(make_pair(prob,GOterm));
        ic_dict[GOterm]=-log(prob);
        
        num_with_parent=allterm_dict.size();
        if (is_a_dict[GOterm].size())
        {
            if (is_a_dict[GOterm].size()==1)
            {
                parent=is_a_dict[GOterm][0];
                num_with_parent=GOterm_dict[parent].size();
            }
            else
            {
                map<string,size_t> accession_count;
                for (i=0;i<is_a_dict[GOterm].size();i++)
                {
                    parent=is_a_dict[GOterm][i];
                    for (a=0;a<GOterm_dict[parent].size();a++)
                    {
                        accession=GOterm_dict[parent][a];
                        if (accession_count.count(accession)==0)
                        {
                            accession_list.push_back(accession);
                            accession_count[accession]=1;
                        }
                        else accession_count[accession]++;
                    }
                }
                num_with_parent=0;
                for (a=0;a<accession_list.size();a++)
                {
                    num_with_parent+=(accession_count[accession_list[a]]
                        ==is_a_dict[GOterm].size());
                    accession_list[a].clear();
                }
                accession_list.clear();
                map<string,size_t> ().swap(accession_count);
            }
        }



        ic_condition_dict[GOterm]=-log(
            (1.+GOterm_dict[GOterm].size())/(1.+num_with_parent));

        if (prob==1) ic_dict[GOterm]=2.22E-16;
        if (GOterm_dict[GOterm].size()==num_with_parent) ic_condition_dict[GOterm]=2.22E-16;
    }

    sort(prob_vec.begin(), prob_vec.end());
    ofstream fout;
    bool toStdout=(outputfilename=="-");
    if (!toStdout) fout.open(outputfilename.c_str());
    for (g=0;g<prob_vec.size();g++)
    {
        prob=prob_vec[prob_vec.size()-g-1].first;
        GOterm=prob_vec[prob_vec.size()-g-1].second;
        if (toStdout) cout<<GOterm<<'\t'<<setiosflags(ios::fixed)
            <<setprecision(6)<<prob<<'\t'<<ic_dict[GOterm]<<'\t'
            <<ic_condition_dict[GOterm]<<'\t'<<name_dict[GOterm]<<endl;
        else fout<<GOterm<<'\t'<<setiosflags(ios::fixed)
            <<setprecision(6)<<prob<<'\t'<<ic_dict[GOterm]<<'\t'
            <<ic_condition_dict[GOterm]<<'\t'<<name_dict[GOterm]<<endl;
    }
    if (!toStdout) fout.close();


    /* clean up */
    map<string,vector<string> > ().swap(allterm_dict);
    vector<string> ().swap(line_vec);
    vector<string> ().swap(accession_list);
    vector<string> ().swap(GOterm_list);
    map<string,vector<string> > ().swap(GOterm_dict);
    vector<pair<double,string> > ().swap(prob_vec);
    map<string,double> ().swap(ic_dict);
    map<string,double> ().swap(ic_condition_dict);
    return;
}

int main(int argc,char **argv)
{
    if (argc!=5)
    {
        cerr<<docstring<<endl;;
        return 0;
    }

    string inputfilename    = argv[1];
    string outputfilename   = argv[2];
    string is_a_filename    = argv[3];
    string name_filename    = argv[4];
    map<string,vector<string> > is_a_dict;
    map<string,string> name_dict;

    parse_is_a_file(is_a_filename, is_a_dict);
    parse_name_file(name_filename, name_dict);
    calculate_ic(inputfilename,outputfilename,is_a_dict,name_dict);

    string ().swap(inputfilename);
    string ().swap(outputfilename);
    string ().swap(is_a_filename);
    string ().swap(name_filename);
    map<string,vector<string> > ().swap(is_a_dict);
    map<string,string> ().swap(name_dict);
    return 0;
}
