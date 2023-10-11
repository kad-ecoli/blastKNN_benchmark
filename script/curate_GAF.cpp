const char* docstring="\n"
"zcat goa_uniprot_all.gaf.gz| curate_GAF - uniprot_sprot_exp [alt_id.csv] [is_a.csv] [evidence.list]\n"
"    generate compact format go annotation\n"
"\n"
"Input:\n"
"    goa_uniprot_all.gaf.gz - GAF format GO annotation in the following format\n"
"                             [1]  DB\n"
"                             [2]  DB Object ID\n"
"                             [3]  DB Object Symbol\n"
"                             [4]  Qualifier, annotations with NOT are ignored\n"
"                             [5]  GO ID\n"
"                             [6]  DB:Reference\n"
"                             [7]  Evidence Code, annotations with ND are ignored\n"
"                             [8]  With (or) From\n"
"                             [9]  Aspect\n"
"                             [10] DB Object Name\n"
"                             [11] DB Object Synonym\n"
"                             [12] DB Object Type\n"
"                             [13] Taxon\n"
"                             [14] Date\n"
"                             [15] Assigned By\n"
"                             [16] Annotation Extension\n"
"                             [17] Gene Product Form ID\n"
"    alt_id.csv             - alt_id GO_ID\n"
"    is_a.csv               - GO_ID Aspect is_a_direct is_a_indirect\n"
"    evidence.list          - evidence code to use\n"
"\n"
"Output:\n"
"    uniprot_sprot_exp.[F,P,C]      - compact go annotation in the format:\n"
"                                     [1] DB Object ID\n"
"                                     [2] Comma separated list of child GO ID\n"
"                                     MF is ignored if the only child term is\n"
"                                     GO:0005515 protein binding\n"
"    uniprot_sprot_exp.list         - list of entries with at least one GO term\n"
"    uniprot_sprot_exp.[F,P,C].is_a - compact go annotation with parent terms\n"
;

#include <iostream>
#include "StringTools.hpp"

using namespace std;

void parse_alt_id_file(const string &alt_id_filename, 
    map<string,string> &alt_id_dict)
{
    bool fromStdin=(alt_id_filename=="-");
    ifstream fin;
    vector<string> line_vec;
    string line;
    if (!fromStdin) fin.open(alt_id_filename.c_str());
    while((fromStdin)?cin.good():fin.good())
    {
        if (fromStdin) getline(cin,line);
        else           getline(fin,line);
        if (line.size()==0) continue;
        Split(line, line_vec, '\t');
        alt_id_dict[line_vec[0]]=line_vec[1];
        line_vec[0].clear(); line_vec[1].clear(); line_vec.clear();
    }
    if (!fromStdin) fin.close();
    vector<string> ().swap(line_vec);
    string ().swap(line);
}

void parse_is_a_file(const string &is_a_filename, 
    map<string,vector<string> >&is_a_dict)
{
    bool fromStdin=(is_a_filename=="-");
    ifstream fin;
    vector<string> line_vec;
    string line;
    string GOterm,Aspect,is_a_direct,is_a_indirect;
    int i;
    if (!fromStdin) fin.open(is_a_filename.c_str());
    while((fromStdin)?cin.good():fin.good())
    {
        if (fromStdin) getline(cin,line);
        else           getline(fin,line);
        if (line.size()==0) continue;
        Split(line, line_vec, '\t');
        GOterm       =line_vec[0];
        Aspect       =line_vec[1];
        is_a_direct  =line_vec[2];
        is_a_indirect=line_vec[3];
        for (i=0;i<line_vec.size();i++) line_vec[i].clear(); line_vec.clear();
        Split(is_a_direct,line_vec,',');
        is_a_dict[GOterm]=line_vec;
        for (i=0;i<line_vec.size();i++) line_vec[i].clear(); line_vec.clear();
        Split(is_a_indirect,line_vec,',');
        for (i=0;i<line_vec.size();i++)
        {
            is_a_dict[GOterm].push_back(line_vec[i]);
            line_vec[i].clear();
        }
        line_vec.clear();
        GOterm.clear();
        Aspect.clear();
        is_a_direct.clear();
        is_a_indirect.clear();
        is_a_dict[GOterm].push_back(GOterm);
    }
    if (!fromStdin) fin.close();
    vector<string> ().swap(line_vec);
    string ().swap(line);
    string ().swap(GOterm);
    string ().swap(Aspect);
    string ().swap(is_a_direct);
    string ().swap(is_a_indirect);
}

void backpropagate(const vector<string> &accession_list, 
    const vector<string>& Aspect_list,
    map<string, map<string, vector<string> > >&GAF_dict,
    map<string,vector<string> > &is_a_dict, const string &prefix)
{
    string txt;
    string Aspect;
    string GOterm,parent;
    string accession;
    vector<string> GOterm_list;
    size_t a,i,j,k;
    ofstream fout;

    for (a=0;a<Aspect_list.size();a++)
    {
        Aspect=Aspect_list[a];
        cout<<prefix+"."+Aspect+".is_a"<<endl;
        txt="";
        for (i=0;i<accession_list.size();i++)
        {
            accession=accession_list[i];
            if (GAF_dict[Aspect].count(accession)==0 ||
                GAF_dict[Aspect][accession].size()==0) continue;
            for (j=0;j<GAF_dict[Aspect][accession].size();j++)
            {
                GOterm=GAF_dict[Aspect][accession][j];
                if (find(GOterm_list.begin(),GOterm_list.end(),GOterm)!=
                     GOterm_list.end()) continue;
                GOterm_list.push_back(GOterm);

                for (k=0;k<is_a_dict[GOterm].size();k++)
                {
                    parent=is_a_dict[GOterm][k];
                    if (find(GOterm_list.begin(),GOterm_list.end(),parent
                        )!=GOterm_list.end()) continue;
                    GOterm_list.push_back(parent);
                }
            }
            sort(GOterm_list.begin(), GOterm_list.end());
            txt+=accession+'\t'+Join(",",GOterm_list)+'\n';
            for (j=0;j<GOterm_list.size();j++) GOterm_list[j].clear();
            GOterm_list.clear();
        }
        fout.open(prefix+"."+Aspect+".is_a");
        fout<<txt<<flush;
        fout.close();
    }
    
    string ().swap(txt);
    vector<string> ().swap(GOterm_list);
    string ().swap(Aspect);
    string ().swap(GOterm);
    string ().swap(parent);
}

void curate_GAF(const string &inputfilename, vector<string>&accession_list,
    vector<string> &Aspect_list, 
    map<string, map<string, vector<string> > > &GAF_dict,
    map<string,string> & alt_id_dict,map<string,vector<string> >&is_a_dict,
    map<string,bool> &evidence_dict, const string &prefix)
{
    bool fromStdin=(inputfilename=="-");
    ifstream fin;
    vector<string> line_vec;
    vector<string> tmp_vec;
    map<string, vector<string> > tmp_dict;
    vector<string> obsolete_list;
    string line;
    string accession,GOterm,evidence,Aspect;
    size_t i;
    if (!fromStdin) fin.open(inputfilename.c_str());
    while((fromStdin)?cin.good():fin.good())
    {
        if (fromStdin) getline(cin,line);
        else           getline(fin,line);
        if (line.size()==0 || line[0]=='!') continue;
        Split(line, line_vec, '\t', false);
        if (line_vec.size()>8 && line_vec[6]!="ND" && 
            line_vec[3].find("NOT") == string::npos)
        {
            accession=line_vec[1];
            GOterm   =line_vec[4];
            evidence =line_vec[6];
            Aspect   =line_vec[8];
            //cout<<accession<<'\t'<<Aspect<<'\t'<<GOterm<<endl;
            if (alt_id_dict.size() && alt_id_dict.count(GOterm))
                GOterm=alt_id_dict[GOterm];
            if (is_a_dict.size() && is_a_dict.count(GOterm)==0)
            {
                if (find(obsolete_list.begin(),obsolete_list.end(),GOterm
                    )==obsolete_list.end())
                {
                    cerr<<"skip "<<GOterm<<endl;
                    obsolete_list.push_back(GOterm);
                }
            }
            else if (evidence_dict.size()    && 
                     evidence_dict.count(evidence)==0); // skip
            else if (evidence_dict.size()==0 && 
                    (evidence=="IBA" || evidence=="ND")); // skip
            else
            {
                if (GAF_dict.count(Aspect)==0)
                {
                    GAF_dict[Aspect]=tmp_dict;
                    Aspect_list.push_back(Aspect);
                }
                if (GAF_dict[Aspect].count(accession)==0)
                    GAF_dict[Aspect][accession]=tmp_vec;
                if (((Aspect=="F" && GAF_dict[Aspect][accession].size()<=3) ||
                     (Aspect!="F" && GAF_dict[Aspect][accession].size()<=1))&&
                    GOterm!="GO:0005515" && GOterm!="GO:0005488" &&
                    GOterm!="GO:0003674" && GOterm!="GO:0008150" &&
                    GOterm!="GO:0005575" && (accession_list.size()==0 ||
                    accession_list.back()!=accession) &&
                    find(accession_list.begin(),accession_list.end(),
                    accession)==accession_list.end())
                    accession_list.push_back(accession);
                if (GAF_dict[Aspect][accession].size()==0 || 
                    find(GAF_dict[Aspect][accession].begin(),
                         GAF_dict[Aspect][accession].end(),GOterm)==
                         GAF_dict[Aspect][accession].end())
                    GAF_dict[Aspect][accession].push_back(GOterm);
            }
        }
        for (i=0;i<line_vec.size();i++) line_vec[i].clear(); line_vec.clear();
    }
    if (!fromStdin) fin.close();
    if (obsolete_list.size()) cerr<<"skip "<<obsolete_list.size()<<" terms"<<endl;
    sort(accession_list.begin(), accession_list.end());

    /* remove entries with only protein binding */
    size_t a;
    for (a=0;a<Aspect_list.size();a++)
    {
        Aspect=Aspect_list[a];
        cout<<GAF_dict[Aspect].size()<<" accession with "<<Aspect<<endl;
    }
    int nonProteinBinding=0;
    if (GAF_dict.count("F"))
    {
        Aspect="F";
        for (a=0;a<accession_list.size();a++)
        {
            accession=accession_list[a];
            if (GAF_dict["F"].count(accession)==0||
                GAF_dict["F"][accession].size()>3) continue;
            nonProteinBinding=0;
            for (i=0;i<GAF_dict["F"][accession].size();i++)
            {
                if (GAF_dict["F"][accession][i]=="GO:0005515" ||
                    GAF_dict["F"][accession][i]=="GO:0005488" ||
                    GAF_dict["F"][accession][i]=="GO:0003674") continue;
                nonProteinBinding=1;
                break;
            }
            if (nonProteinBinding==0) GAF_dict[Aspect].erase(accession);
        }
        cout<<GAF_dict[Aspect].size()<<" accession with "<<Aspect<<endl;
    }

    /* output */
    ofstream fout;
    string txt;
    for (i=0;i<accession_list.size();i++)
        txt+=accession_list[i]+"\n";
    cout<<prefix+".list"<<endl;
    fout.open(prefix+".list");
    fout<<txt<<flush;
    fout.close();
    for (a=0;a<Aspect_list.size();a++)
    {
        if (is_a_dict.size()) break;
        Aspect=Aspect_list[a];
        cout<<prefix+"."+Aspect<<endl;
        txt="";
        for (i=0;i<accession_list.size();i++)
        {
            accession=accession_list[i];
            if (GAF_dict[Aspect].count(accession) && 
                GAF_dict[Aspect][accession].size())
            {
                sort(GAF_dict[Aspect][accession].begin(),
                     GAF_dict[Aspect][accession].end());
                txt+=accession+'\t'+Join(",",GAF_dict[Aspect][accession])+'\n';
            }
        }
        if (txt.size())
        {
            fout.open(prefix+"."+Aspect);
            fout<<txt<<flush;
            fout.close();
        }
    }


    /* clean up */
    vector<string> ().swap(line_vec);
    string ().swap(line);
    string ().swap(txt);
    map<string, vector<string> > ().swap(tmp_dict);
    vector<string> ().swap(tmp_vec);
    vector<string> ().swap(obsolete_list);
}

void parse_evidence_file(const string &evidence_filename, 
    map<string,bool> &evidence_dict)
{
    bool fromStdin=(evidence_filename=="-");
    ifstream fin;
    string line;
    if (!fromStdin) fin.open(evidence_filename.c_str());
    while((fromStdin)?cin.good():fin.good())
    {
        if (fromStdin) getline(cin,line);
        else           getline(fin,line);
        if (line.size()==0) continue;
        evidence_dict[line]=1;
    }
    if (!fromStdin) fin.close();
    string ().swap(line);
}

int main(int argc,char **argv)
{
    if (argc<3 || argc>6)
    {
        cerr<<docstring<<endl;;
        return 0;
    }

    string inputfilename    = argv[1];
    string prefix           = argv[2];
    string alt_id_filename  =(argc>3)?argv[3]:"";
    string is_a_filename    =(argc>4)?argv[4]:"";
    string evidence_filename=(argc>5)?argv[5]:"";
    map<string,string> alt_id_dict;
    map<string,vector<string> > is_a_dict;
    map<string,bool> evidence_dict;

    if (alt_id_filename.size())
    {
        parse_alt_id_file(alt_id_filename, alt_id_dict);
        if (is_a_filename.size())
        {
            parse_is_a_file(is_a_filename, is_a_dict);
            if (evidence_filename.size())
                parse_evidence_file(evidence_filename, evidence_dict);
        }
    }
    vector<string> accession_list;
    vector<string> Aspect_list;
    map<string, map<string, vector<string> > > GAF_dict;
    map<string, vector<string> > tmp_dict;
    curate_GAF(inputfilename,accession_list,Aspect_list,
        GAF_dict,alt_id_dict,is_a_dict,evidence_dict,prefix);
    if (is_a_filename.size())
        backpropagate(accession_list,Aspect_list,GAF_dict,is_a_dict,prefix);

    string ().swap(inputfilename);
    string ().swap(is_a_filename);
    string ().swap(alt_id_filename);
    string ().swap(prefix);
    map<string,string> ().swap(alt_id_dict);
    map<string,vector<string> > ().swap(is_a_dict);
    vector<string> ().swap(accession_list);
    vector<string> ().swap(Aspect_list);
    map<string, map<string, vector<string> > >().swap(GAF_dict);
    map<string,bool> ().swap(evidence_dict);
    return 0;
}
