#include "writeme.h"

writeme::writeme(string mainpath)
{
    main_path = mainpath;
}

writeme::writeme(string mainpath, string default_filename)
{
    main_path = mainpath;
    filename = default_filename;
    using_default_filename = true;
}

writeme::~writeme()
{
}

void writeme::set_subpath(string subpath)
{
    sub_path = subpath;
    using_sub_directory = true;
}

void writeme::clear_subpath()
{
    sub_path.clear();
    using_sub_directory=false;
}

void writeme::set_default_filename(string filename_, bool add_ints, bool add_doubles)
{
    filename.clear();
    filename = filename_;
    if (add_ints || add_doubles)
    {
        filename += parameter_string(add_ints,add_doubles);
        using_values_at_end = true;
    }
    using_default_filename = true;
}

void writeme::clear_default_filename()
{
    filename.clear();
    using_default_filename = false;
    using_values_at_end = false;
}

void writeme::add_double(string key, double value)
{
    doubles[key] = value;
}

void writeme::remove_double(string key)
{
    doubles.erase(key);
}

void writeme::add_int(string key, int value)
{
    ints[key] = value;
}

void writeme::remove_int(string key)
{
    ints.erase(key);
}

string writeme::to_string_no_trailing(double value)
{
    string temp;
    temp = to_string(value);
    if (temp.find('.') == string::npos)
    {
        return temp;
    }
    while (temp.back() == '0')
    {
        temp.pop_back();
        if (temp.back() == '.')
        {
            temp.pop_back();
            break;
        }
    }
    return temp;
}

string writeme::to_string_no_trailing(int value)
{
    return to_string_no_trailing(double(value));
}

string writeme::parameter_string(bool add_ints, bool add_doubles)
{
    string temp;
    if (add_doubles)
    {
        for (auto &key : doubles)
        {
            temp += "_" + key.first + "_" + to_string_no_trailing(key.second);
        }
    }

    if (add_ints)
    {
        for (auto &key : ints)
        {
            temp += "_" + key.first + "_" + to_string_no_trailing(key.second);
        }
    }
    return temp;
}

string writeme::get_current_path(string filenamestart, bool add_ints, bool add_doubles)
{
    string temp = main_path;
    if (using_sub_directory)
    {
        temp += sub_path;
    }

    if (filenamestart != "")
    {
        temp += filenamestart;
    }
    if (using_default_filename)
    {
        temp += filename;
    }
    if (!using_values_at_end && (add_ints || add_doubles))
    {

        temp += parameter_string(add_ints,add_doubles);
    }
    return temp;
}

void writeme::print_current_path(string filenamestart)
{
    cout << get_current_path(filenamestart) << endl;
}


//------------------------------------------------------------------------------------
//                                  WRITE FUNCTIONS
//------------------------------------------------------------------------------------

void writeme::write_double_vector(vector<double> quantity,string filenamestart, bool add_ints, bool add_doubles){

    string temp;
    temp = get_current_path(filenamestart,add_ints,add_doubles) + ".txt";

    cout << "Writing \n";
    cout << temp << endl;
    ofstream file(temp);
    for (const auto &e : quantity) file << setw(15) << setprecision(8) << e << "\n";
    file.close();
}

void writeme::write_vector_vector(vector<vector<double>> quantity, string filenamestart, bool add_ints, bool add_doubles)
{
    string temp;
    temp = get_current_path(filenamestart, add_ints, add_doubles) + ".txt";

    cout << "Writing \n";
    cout << temp << endl;
    ofstream file(temp);
    file << setiosflags(ios::showpoint | ios::uppercase);

    for (u_long j = 0; j < quantity[0].size(); ++j){
        for (u_long i = 0; i < quantity.size(); ++i) {
            file << setw(15) << setprecision(10) << quantity[i][j];
        }
        file << "\n";
    }
    file.close();
}

void writeme::write_vector_vector_bin(vector<vector<double>> quantity, string filenamestart, bool add_ints, bool add_doubles){
    string temp;
    temp = get_current_path(filenamestart, add_ints, add_doubles) + ".bin";
    cout << "Writing \n";
    cout << temp << endl;
    ofstream file(temp, ofstream::binary);
    for (u_long i = 0; i < quantity.size(); ++i) {
        file.write(reinterpret_cast<char*>(&quantity[i][0]),quantity[i].size()*sizeof(double));

    }
    file.close();
}
