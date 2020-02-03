#ifndef WRITEME_H
#define WRITEME_H

#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>

/*
 *
 * WORK IN PROGRESS
 * this class is inteded for further development for printing of many types



TODO:
Add option to add parameters in first line of file for textfiles
check for existing path
Add workexample


*/

using namespace std;

class writeme
{
private:
    // Path stuff
    string main_path;

    bool using_sub_directory = false;
    string sub_path;

    // Filenamestuff
    string filename;
    map<string,double> doubles;
    map<string,int> ints;
    bool using_default_filename = false;
    bool using_values_at_end = false;


public:
    writeme(string mainpath);
    writeme(string mainpath, string default_filename);
    ~writeme();

    // Config
    string parameter_string(bool add_ints = false, bool add_doubles = false);
    void set_default_filename(string filename_,bool add_ints = false, bool add_doubles = false);
    void clear_default_filename();
    void set_subpath(string subpath);
    void clear_subpath();
    void add_double(string key, double value);
    void remove_double(string key);
    void add_int(string key, int value);
    void remove_int(string key);
    void clear_maps() {doubles.clear(); ints.clear();}

    // Convenient
    string get_current_path(string filenamestart = "", bool add_ints = false, bool add_doubles = false);
    void print_current_path(string filenamestart = "");
    string to_string_no_trailing(double value);
    string to_string_no_trailing(int value);

    // Writing
    void write_double_vector(vector<double> quantity,string filenamestart = "",bool add_ints = false, bool add_doubles = false);
    void write_vector_vector(vector<vector<double>> quantity,string filenamestart = "", bool add_ints = false, bool add_doubles = false);
    void write_vector_vector_bin(vector<vector<double>> quantity, string filenamestart = "", bool add_ints = false, bool add_doubles = false);
};

#endif // WRITEME_H
