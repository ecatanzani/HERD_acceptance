#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <sstream>

#include "TTree.h"
#include "TFile.h"

#define words_in_line 54
#define glat_word_counter 23

std::istringstream readStreamFromFile(std::string input_RTI_path)
{
    std::ifstream input_file(input_RTI_path.c_str());
    if (!input_file.is_open())
    {
        std::cerr << "\nERROR 100! File not open " << input_RTI_path << "\n\n";
        exit(100);
    }
    std::string input_string((std::istreambuf_iterator<char>(input_file)), (std::istreambuf_iterator<char>()));
    input_file.close();
    std::istringstream input_stream(input_string);
    return input_stream;
}

void parse_stream(
    const std::string input_RTI_path,
    double &glat, 
    double &glon, 
    TTree &RTI_tree)
{
    auto input_stream = readStreamFromFile(input_RTI_path);
    std::string tmp_RTI_file;

    while (input_stream >> tmp_RTI_file)
    {
        auto rti_stream = readStreamFromFile(tmp_RTI_file);
        // Create line counter variable
        unsigned int line_counter = 0;

        std::string tmp_str;
        std::string header;
        std::string::size_type sz;

        // parse stream
        //while (input_stream >> tmp_str)
        rti_stream >> header;
        
        while(!rti_stream.eof())
        {
            rti_stream >> tmp_str;
            auto JMDCtime = stoul(tmp_str, &sz);
            rti_stream >> tmp_str;
            auto run = stoul(tmp_str, &sz);
            if (!run)
            {
                ++line_counter;
                continue;
            }

            // skipping first glat_word_counter -1 words and read the glat_word_counter-th word
            for (auto wIdx=0; wIdx<(glat_word_counter-2); ++wIdx)
                rti_stream >> tmp_str;
            
            glat = stod(tmp_str, &sz);
            rti_stream >> tmp_str;
            glon = stod(tmp_str, &sz);
            
            std::cout << "\nLat: " << glat << "\tLon: " << glon << std::endl;

            // Finish reading RTI line
            for (auto wIdx=0; wIdx<(words_in_line-(glat_word_counter+1)); ++wIdx)
                rti_stream >> tmp_str;
            
            ++line_counter;
            
            RTI_tree.Fill();
        }
        
        std::cout << "\nHas been read " << lines << " lines from: " << tmp_RTI_file << std::endl;
    }
}

void csvToTree(std::string input_RTI_path, std::string outFile)
{
    // Create glon e glat variables
    double glat = -1;
    double glon = -1;

    // Create new TFile
    TFile RTI_output_file(outFile.c_str(),"RECREATE");
    if (RTI_output_file.IsZombie())
    {
        std::cerr << "\n\nError writing output TFile: " << outFile << std::endl;
        exit(100);
    }

    // Create final Tree
    TTree RTI_tree("RTI_tree", "AMS RTI TTree");
    // Branch the TTree
    RTI_tree.Branch("glat", &glat, "glat/D");
    RTI_tree.Branch("glon", &glon, "glon/D");

    // Parse stream
    parse_stream(
        input_RTI_path,
        glat, 
        glon, 
        RTI_tree);
    
    RTI_tree.Write();
}