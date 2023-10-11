const char* docstring="\n"
"subset_tsv input.is_a input.list output.F.is_a\n"
"    extract a subset of GOA entries\n"
"\n"
"Input:\n"
"    input.is_a - COFACTOR format GO annotation.\n"
"                 [1] accession\n"
"                 [2] comma separated list of GO terms\n"
"    input.list - list of accession to parse\n"
"\n"
"Output:\n"
"    output.is_a - subset of GO annotation listed by input.list\n"
;

#include <iostream>
#include "StringTools.hpp"

using namespace std;

void read_accession_list(const string &inputListFile,
    vector<string> &accession_list)
{
    bool fromStdin=(inputListFile=="-");
    ifstream fin;
    string line;
    if (!fromStdin) fin.open(inputListFile.c_str());
    while((fromStdin)?cin.good():fin.good())
    {
        if (fromStdin) getline(cin,line);
        else           getline(fin,line);
        if (line.size()==0) continue;
        accession_list.push_back(line);
        line.clear();
    }
    if (!fromStdin) fin.close();
}

void subset_tsv(const string &inputfilename, 
    const vector<string>&accession_list, const string &outputfilename)
{
    map<string,bool> accession_map;
    size_t i;
    for (i=0;i<accession_list.size();i++)
        accession_map[accession_list[i]]=false;

    bool fromStdin=(inputfilename=="-");
    bool toStdout=(outputfilename=="-");
    ifstream fin;
    ofstream fout;
    vector<string> line_vec;
    string line;
    if (!fromStdin) fin.open(inputfilename.c_str());
    if (!toStdout) fout.open(outputfilename.c_str());
    while((fromStdin)?cin.good():fin.good())
    {
        if (fromStdin) getline(cin,line);
        else           getline(fin,line);
        if (line.size()==0 || line[0]=='!') continue;
        Split(line, line_vec, '\t', false);
        if (line_vec.size()>=2 && 
            accession_map.count(line_vec[0]))
        {
            if (toStdout) cout<<line<<endl;
            else          fout<<line<<endl; 
        }
        for (i=0;i<line_vec.size();i++) line_vec[i].clear();
        line_vec.clear();
    }
    if (!fromStdin) fin.close();
    if (!toStdout) fout.close();

    /* clean up */
    string ().swap(line);
    vector<string>  ().swap(line_vec);
    map<string,bool>().swap(accession_map);
}

int main(int argc,char **argv)
{
    if (argc<4)
    {
        cerr<<docstring<<endl;;
        return 0;
    }

    string inputfilename  = argv[1];
    string inputListFile  = argv[2];
    string outputfilename = argv[3];

    vector<string> accession_list;
    read_accession_list(inputListFile, accession_list);
    subset_tsv(inputfilename, accession_list, outputfilename);

    string ().swap(inputfilename);
    string ().swap(inputListFile);
    string ().swap(outputfilename);
    vector<string> ().swap(accession_list);
    return 0;
}
