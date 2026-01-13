#pragma once

#include<iostream>
#include<iomanip>
#include<string>
#include<array>
#include<vector>
#include<memory>
#include<stdexcept>
#include<sys/stat.h>
#include<Standard_Algorithms/Print_Routines.h>
#include"FM_Error_Handling.h"

namespace File_Management
{

namespace print = Print_Routines;
namespace error = File_Management::Error_Handling;

// function : returns the path to the root folder of the project, which is marked by the .orientation file
std::string Orientation()
{
    std::string tree_climb = "";
    struct stat buffer;
    while (stat((tree_climb + ".orientation").c_str(), &buffer) != 0)
    {
        tree_climb += "../";
        if (tree_climb.length() > 100) // safety check to avoid infinite loop
        {
            throw std::runtime_error("Orientation file not found within reasonable directory depth.");
        }
    }
    return tree_climb;
}

// function : creates a single folder and returns the name (without / in the end)
std::string create_folder( std::string folder_name, const size_t max_name_length = 200 )
{
    print::cut_if_too_large( folder_name, max_name_length );
    char const *c = folder_name.data();
    if( mkdir( c, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH ) != -1 ) // then folder can be built
    { 
        std::cout << "created folder : " << folder_name << "\n";
    }
    return folder_name;
}

// function : creates a folder tree and returns the total name (without / in the end)
std::string create_folder_tree( std::vector<std::string> folders, const size_t max_name_length = 200 )
{
    std::string folder_name = "";
    for( auto & subfolder : folders )
    {
        folder_name += subfolder;
        create_folder( folder_name, max_name_length );
        folder_name += "/";
    }
    return folder_name;
}

// obtain the system name by executing uname -n
std::string get_system_name()
{
    const char* cmd{ "uname -n" };
    std::array<char, 128> buffer;
    std::string result{};
    std::unique_ptr<FILE, int(*)(FILE*)> pipe(popen(cmd, "r"), pclose);
    if( !pipe )
    {
        throw std::runtime_error("popen() failed!");
    }
    while( fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr )
    {
        result += buffer.data();
    }

    // remove the trailing newline character
    if( !result.empty() && result[result.size() - 1] == '\n' )
    {
        result.erase(result.size() - 1);
    }
    return result;
}

};