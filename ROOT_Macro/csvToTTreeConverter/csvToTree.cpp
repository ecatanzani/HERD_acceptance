#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <sstream>

#include "TTree.h"
#include "TFile.h"
#include "TMath.h"

#define words_in_line 54
#define glat_word_counter 23
#define thetas_word_counter 19

std::string readStreamFromFile(std::string input_RTI_path)
{
    std::ifstream input_file(input_RTI_path.c_str());
    if (!input_file.is_open())
    {
        std::cerr << "\nERROR 100! File not open " << input_RTI_path << "\n\n";
        exit(100);
    }
    std::string input_string((std::istreambuf_iterator<char>(input_file)), (std::istreambuf_iterator<char>()));
    input_file.close();
    return input_string;
}

void parse_stream(
    const std::string input_RTI_path,
    double &glat, 
    double &glon,
    double &thetas,
    double &phi,
    int &events,
    double &geo_lat,
    double &geo_lon,
    TTree &RTI_tree)
{
    std::istringstream input_stream(readStreamFromFile(input_RTI_path));
    std::string tmp_RTI_file;

    while (input_stream >> tmp_RTI_file)
    {
        std::istringstream rti_stream(readStreamFromFile(tmp_RTI_file));
        // Create line counter variable
        unsigned int line_counter = 0;
        
        std::string tmp_str;
        std::string header;
        std::string::size_type sz;

        // parse stream
        //while (input_stream >> tmp_str)
        rti_stream >> header;
        int wcounter = 0;

        while(rti_stream >> tmp_str)
        {
            ++wcounter;

            if (wcounter == 1)
                auto JMDCtime = stoul(tmp_str, &sz);
            else if (wcounter == 2)
            {
                auto run = stoul(tmp_str, &sz);
                if (!run)
                {
                    ++line_counter;
                    wcounter=0;
                    continue;
                }
            }
            else if (wcounter == thetas_word_counter)
                thetas = stod(tmp_str, &sz);
            else if (wcounter == (thetas_word_counter+1))
                phi = stod(tmp_str, &sz);
            else if (wcounter == glat_word_counter) 
                glat = stod(tmp_str, &sz);
            else if (wcounter == (glat_word_counter+1))
                glon = stod(tmp_str, &sz);
            else if (wcounter == (glat_word_counter+2))
                events = stoi(tmp_str, &sz);
            else if (wcounter == words_in_line)
            {
                wcounter=0;
                ++line_counter;
                geo_lat = (1-thetas)*TMath::RadToDeg();
                geo_lon = phi*TMath::RadToDeg();
                RTI_tree.Fill();
            }
        }
        
        std::cout << "\nHas been read " << line_counter << " lines from: " << tmp_RTI_file;
    }
}

void csvToTree(std::string input_RTI_path, std::string outFile)
{
    // Create glon e glat variables
    double glat = -1;
    double glon = -1;
    double thetas = -1;
    double phi = -1;
    int events = -1;
    double geo_lat = -1;
    double geo_lon = -1;

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
    RTI_tree.Branch("thetas", &thetas, "thetas/D");
    RTI_tree.Branch("phi", &phi, "phi/D");
    RTI_tree.Branch("events", &events, "events/I");
    RTI_tree.Branch("geo_lat", &geo_lat, "geo_lat/D");
    RTI_tree.Branch("geo_lon", &geo_lon, "geo_lon/D");

    // Parse stream
    parse_stream(
        input_RTI_path,
        glat, 
        glon,
        thetas,
        phi,
        events,
        geo_lat,
        geo_lon,
        RTI_tree);
    
    RTI_tree.Write();
}