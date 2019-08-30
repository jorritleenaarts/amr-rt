#include <fstream>      // std::ifstream
#include <string>
#include <fstream>
#include <sstream>
#include <dirent.h>
#include <sys/stat.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cstring>
#include <regex>
#include <typeinfo>

#include "ftools.hpp"

using namespace std;

/* ---
   
   File and string manipulation tools.
   Also my own parser for input files: A hybrid between M. van Noort 
   implementation and the usual [Section] with values.

   Some string manipulation routines adapted from Stackoverflow.

   Coded by J. de la Cruz Rodriguez (ISP-SU 2017/2018)
   
   ---- */


/* ---------------------------------------------------- */

bool ftl::file_exists(const string name) {
    ifstream f(name.c_str(),ifstream::in);
    if (f.good()) {
        f.close();
        return true;
    } else {
        f.close();
        return false;
    }
}

/* ---------------------------------------------------- */

vector<string> ftl::readFile(string filein, const char com, bool remove_empty)
{

    bool exists = file_exists(filein);
    if(!exists) fprintf(stderr,"error: ftl::readFile: cannot open file [%s]\n", filein.c_str());

    ifstream in(filein.c_str(), ios::in | ios::binary);
    vector<string> res;
    string line, com1(" ");
    com1.assign(1, com);

    if(in){
        while(getline(in, line)){
            if(remove_empty) if(line == "") continue;
            if(line[0] != com) res.push_back(ftl::cleanLine(line,com1));
        }
    }

    return res;

}

// *********************************************************************

ftl::section_t ftl::readConfigFile(const string& config_file)
{
    vector<string> txt = ftl::readFile(config_file);
    size_t off = 0;
    ftl::section_t configuration = ftl::parseSection(txt, off, 1);
    return configuration;
}

/* ---------------------------------------------------- */

string ftl::toUpper(string in)
{
    transform(in.begin(), in.end(), in.begin(), ::toupper);
    return in;
}

/* ---------------------------------------------------- */

string ftl::toLower(string in)
{
    transform(in.begin(), in.end(), in.begin(), ::tolower);
    return in;
}

/* ---------------------------------------------------------------------------- */

string ftl::removeSpaces(string input){
    input.erase(remove(input.begin(),input.end(),' '),input.end());
    return input;
}

/* ---------------------------------------------------------------------------- */

vector<string> ftl::strsplit( const string &var, string token, bool rmspaces){


    vector<string> res;
    stringstream ss(var);

    string tmp;
    if(rmspaces) while(getline(ss, tmp, *(token.c_str()))){
            tmp = removeSpaces(tmp);
            if(tmp.length() > 0) res.push_back(tmp);
        }else
        while(getline(ss, tmp, *(token.c_str()))){
            res.push_back(tmp);
        }

    return res;
}

/* ---------------------------------------------------------------------------- */

vector<string> ftl::strsplit2(const string& s, string rgx_str = "\\s+") {

    vector<string> elems;
    regex rgx (rgx_str);

    sregex_token_iterator iter(s.begin(), s.end(), rgx, -1);
    sregex_token_iterator end;

    while (iter != end)  {
        if (iter->length() > 0) elems.push_back(*iter);
        ++iter;
    }

    return elems;

}

/* ---------------------------------------------------------------------------- */

bool ftl::bdir_exists(string name){
    DIR* dir = opendir(name.c_str());
    if (dir) return true;
    else if (ENOENT == errno) return false;
    else{
        fprintf(stderr,"error: bdir_exists: cannot write to folder [%s], exiting.\n", name.c_str());
        return true; // Just ensure function will return a value
    }
}

/* ---------------------------------------------------------------------------- */

string ftl::cleanLine(string var,  string token, bool rmspaces)
{
    vector<string> res = ftl::strsplit(var, token, rmspaces);
    return res[0];
}

/* ---------------------------------------------------------------------------- */

string ftl::trimString(const string str, const string whitespace)
{
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

/* ---------------------------------------------------------------------------- */


string ftl::reduceString(const string str, const string fill,
                         const string whitespace)
{
    // trim first
    auto result = ftl::trimString(str, whitespace);

    // replace sub ranges
    auto beginSpace = result.find_first_of(whitespace);
    while (beginSpace != string::npos)
    {
        const auto endSpace = result.find_first_not_of(whitespace, beginSpace);
        const auto range = endSpace - beginSpace;

        result.replace(beginSpace, range, fill);

        const auto newStart = beginSpace + fill.length();
        beginSpace = result.find_first_of(whitespace, newStart);
    }

    return result;
}

/* ---------------------------------------------------------------------------- */

ftl::section_t ftl::parseSection(vector<string> &txt, size_t &off_in,  int is_main)
{
    const string tok = "=";

    section_t sec;
    ftl::field_t field;
    if(is_main > 0) sec.ID = "[/]";


    int in_scope = 0;
    size_t off = (size_t)off_in, nlines = (size_t)txt.size();

    while(off < nlines){

        string tmp = ftl::trimString(ftl::cleanLine(txt[off], "#", false));
        unsigned tmp_len = (unsigned)tmp.size();
        if((tmp[0] == '[') && (tmp[tmp_len-1] == '{')){
            if(in_scope == 0 && (is_main == 0)){
                in_scope += 1;
                sec.ID = tmp.substr(0,tmp.size()-1);
                sec.fields.clear();
            }else sec.sub.push_back(parseSection(txt, off, 0));
        }else if((tmp[0] == '}' && (is_main == 0))){
            off_in = off;
            in_scope -= 1;
            return sec;
        }else{
            vector<string> spl = ftl::strsplit(tmp, tok, false);
            size_t fsiz = spl.size();

            if(fsiz != 2){
                if(fsiz == 1)
                    fprintf(stderr,"warning: ftl::parseArr: ignoring line [%ld] -> [%s]\n", off, tmp.c_str());
                off++;
                continue;
            }

            /* ---- remove all spaces from string --- */

            spl[0] = ftl::reduceString(spl[0]);
            spl[1] = ftl::reduceString(spl[1]);


            field.name = spl[0], field.value = spl[1];
            sec.fields.push_back(field);
        }

        off++;
    }// off

    off_in = off;
    return sec;

}

/* ---------------------------------------------------------------------------- */

bool ftl::searchSection(const ftl::section_t &sec, string sname, ftl::section_t &res, int search_sub)
{

    bool found = false;

    if(sec.ID == sname){
        found = true;
        res = sec;
        return found;
    }else{

        const vector<section_t> &se = sec.sub;
        size_t nsec = (size_t)se.size();

        for(size_t ii=0; ii<nsec; ii++){

            if(se[ii].ID == sname){
                res.ID = se[ii].ID;
                res.fields = se[ii].fields;
                res.sub    = se[ii].sub;
                found = true;
                break;
            }

            if((search_sub > 0) && (se[ii].sub.size() > 0)){
                found = ftl::searchSection(se[ii], sname, res, (int)(search_sub-1));
                if(found) break;
            }

        }// ii
    }
    return found;

}

/* ---------------------------------------------------------------------------- */


ftl::section_t ftl::parseArr(vector<string> &txt, bool verbose)
{
    size_t off = 0;
    section_t sec = ftl::parseSection(txt, off, 1);

    return sec;
}

/* ---------------------------------------------------------------------------- */

bool ftl::searchField(const ftl::section_t &sec, string fname, ftl::field_t &res)
{

    size_t nfields = sec.fields.size();
    bool found = false;

    for(size_t ii=0;ii<nfields;ii++){
        if(sec.fields[ii].name == fname){
            found = true;
            res.name = sec.fields[ii].name;
            res.value = sec.fields[ii].value;
            break;
        }
    }

    return found;
}

void ftl::message(const string messg, const int error_level) {
    if (error_level == 0){
        cout << "MESSAGE: " << messg << endl;
    } else if (error_level == 1) {
        cout << "WARNING: " << messg << endl;
    } else if (error_level == 2) {
        cout << "ERROR: " << messg << endl;
        exit(1);
    }
}

/* ---------------------------------------------------------------------------- */

/*
void ftl::printSection(const ftl::section_t &sec, lg::log &l, size_t ilev)
{

  static const lg::msg_t mtyp = lg::MSGXTRA;
  
  char buf[250];
  size_t nfields = (size_t)sec.fields.size();
  size_t nsub    = (size_t)sec.sub.size();
  
  size_t nilev   = ilev * 3;
  size_t nspace  = nilev + 3;
  
  char space[nspace+1] = {}; memset(space,' ', nspace), space[nspace] = '\0';
  char spacl[ nilev+1] = {}; memset(spacl,' ', nilev),  spacl[nilev] = '\0';
  
  
  
  // --- print this section and subfieds --- //
  
  sprintf(buf, "%s%s{", spacl, sec.ID.c_str());
  l.msg(buf, mtyp);
  
  
  for(size_t ii = 0; ii < nfields; ii++){
    const ftl::field_t &f =  sec.fields[ii];
    sprintf(buf, "%s%s = %s", space, f.name.c_str(), f.value.c_str());
    l.msg(buf, mtyp);
  }



  // --- If there are more sub-sections, print them recursively --- //

  if(nsub > 0)
    for(size_t ii = 0; ii < nsub; ii++)
      ftl::printSection(sec.sub[ii], l, ilev+1);

  sprintf(buf, "%s}", spacl);
  l.msg(buf, mtyp);
}
*/


/* ---------------------------------------------------------------------------- */
