#ifndef __FTL_H__
#define __FTL_H__

/* ---
   
   File and string manipulation tools.
   Also my own parser for input files: A hybrid between M. van Noort 
   implementation and the usual [Section] with values.

   Some string manipulation routines adapted from Stackoverflow.

   Coded by J. de la Cruz Rodriguez (ISP-SU 2017/2018)
   
   
   ---- */

#include <string>
#include <vector>
#include <typeinfo>
#include <algorithm>
#include <iostream>
#include <sstream>

namespace ftl {

    // ---- Types ---- //

    struct field_t {
        std::string name;
        std::string value;
    };

    struct section_t {
        std::string ID;
        std::vector<field_t> fields;
        std::vector<section_t> sub;
    };


    // ---- Prototypes ---- //

    bool file_exists(const std::string name);

    std::vector<std::string> readFile(std::string filein, const char com = '#', bool remove_empty = true);

    section_t readConfigFile(const std::string& config_file);

    std::string toUpper(std::string in);

    std::string toLower(std::string in);

    std::string removeSpaces(std::string input);

    std::vector<std::string> strsplit(const std::string &var, std::string token, bool rmspaces);

    std::vector<std::string> strsplit2(const std::string &s, std::string rgx_str);

    bool bdir_exists(std::string name);

    std::string cleanLine(std::string var, std::string token, bool rmspaces = false);

    section_t parseArr(std::vector<std::string> &txt, bool verbose = false);

    std::string reduceString(const std::string str,
                             const std::string fill = " ",
                             const std::string whitespace = " \t");

    std::string trimString(const std::string str,
                           const std::string whitespace = " \t");

    section_t parseSection(std::vector<std::string> &txt, size_t &off_in, int in_main);

    bool searchSection(const section_t &sec, std::string sname, section_t &res, int search_sub);

    bool searchField(const section_t &sec, std::string fname, field_t &res);

    void message(const std::string messg, const int error_level);

    //  void printSection(const section_t &sec, lg::log &l, size_t ilev = 0);


    // --- Implementation / Templates --- //


    // ---------------------------------------------------- //

    template<typename T>
    bool getKeyword(const ftl::section_t &inputsec,
                   const std::string &secname, const std::string &fieldname, T& var) {

        bool found = true;

        ftl::section_t sec;
        if (!ftl::searchSection(inputsec, secname, sec, 1)) {
            ftl::message("Section named '" + secname + "' does not exist.", 1);
            found = false;
        }

        ftl::field_t res;
        if (found) {
            if (!ftl::searchField(sec, fieldname, res)) {
                std::cout << "Section '" + secname + "' does not contain field '"
                             + fieldname + "'." << std::endl;
                found = false;
            }
        }

        if (!found) {
            std::cout << "Falling back on default value for " + secname + ":"
                         + fieldname+": " << var << "." << std::endl;
        } else {
            std::stringstream ss(res.value);
            ss >> var;
            if (ss.fail()) {
                std::cout << "Error reading keyword " + secname + ":"+fieldname
                             + ". Stopping execution." << std::endl;
                exit(1);
            }
        }

        return found;

    }


    template <class T> T parseValue(field_t &f)
    {
        T val = {};
        char form[7];

        if     (typeid(T) == typeid(double))    sprintf(form, "%s",  "%lf");
        else if(typeid(T) == typeid(float))     sprintf(form, "%s",  "%f");
        else if(typeid(T) == typeid(int))       sprintf(form, "%s",  "%d");
        else if(typeid(T) == typeid(long int))  sprintf(form, "%s",  "%ld");
        else if(typeid(T) == typeid(unsigned))  sprintf(form, "%s",  "%u");
        else if(typeid(T) == typeid(size_t))    sprintf(form, "%s",  "%lu");
        else{
            std::string tmp = typeid(T).name();
            fprintf(stderr,"Warning: ftl::parseValue: var type [%s] not implemented -> [%s] = 0, instead of [%s]\n", tmp.c_str(), f.name.c_str(), f.value.c_str());

            return val;
        }


        // --- read value --- //

        sscanf(f.value.c_str(), form, &val);

        return val;
    }

    // ---------------------------------------------------- //

    template <class T> bool parseField(section_t &sec, std::string section,
                                       std::string field, T &val)
    {

        // --- First check whether the field is in the root struct --- //

        bool found = false;
        //v = 0.0;

        if(sec.ID == section){
            size_t nfields = (size_t)sec.fields.size();

            for(size_t ii=0; ii<nfields; ii++)
                if(sec.fields[ii].name == field){
                    val = parseValue<T>(sec.fields[ii]);
                    return true;
                }

        }else{

            size_t nsub = sec.sub.size();
            for(size_t ii = 0; ii<nsub; ii++){
                if((found = parseField<T>(sec.sub[ii], section, field, val))) return true;
            }

        }


        return found;
    }
};


#endif
