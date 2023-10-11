const char* docstring=""
"tsv2fasta metaclust.pfam metaclust.fasta\n"
"    convert tab-eliminated table metaclust.pfam\n"
"    to FASTA format alignment metaclust.fasta\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>

using namespace std;

int tsv2fasta(const string infile="-", const string outfile="-")
{
    ifstream fp_in;
    ofstream fp_out;
    if (infile!="-") fp_in.open(infile.c_str(),ios::in);
    if (outfile!="-") fp_out.open(outfile.c_str(),ofstream::out);
    string txt,line;
    int nseqs=0;
    int i,j;
    while ((infile!="-")?fp_in.good():cin.good())
    {
        if (infile!="-") getline(fp_in,line);
        else getline(cin,line);

        for (i=0;i<line.size();i++)
        {
            if (line[i]=='\t')
            {
                nseqs++;
                txt='>'+line.substr(0,i)+'\n';
                line=line.substr(i+1);
                txt+=line;
                if (outfile!="-") fp_out<<txt<<endl;
                else                cout<<txt<<endl;
                txt.clear();
                break;
            }
        }
    }
    fp_in.close();
    fp_out.close();
    return nseqs;
}

int main(int argc, char **argv)
{
    /* parse commad line argument */
    if(argc<2)
    {
        cerr<<docstring;
        return 0;
    }
    string infile=argv[1];
    string outfile=(argc<=2)?"-":argv[2];
    tsv2fasta(infile,outfile);
    return 0;
}
