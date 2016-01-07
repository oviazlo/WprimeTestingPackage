#ifndef HelpFunctions_H
#define HelpFunctions_H

#include "TStyle.h"
#include "TROOT.h"
#include <iostream>


int parseOptionsWithBoost(po::variables_map &vm, po::options_description desc, 
                          int argc, char* argv[]);
TStyle* AtlasStyle();
void SetAtlasStyle();

#endif