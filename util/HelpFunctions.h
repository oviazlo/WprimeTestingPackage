#ifndef HelpFunctions_H
#define HelpFunctions_H

#include "TStyle.h"
#include "TROOT.h"
#include <iostream>

namespace 
{ 
  const size_t ERROR_IN_COMMAND_LINE = 1; 
  const size_t SUCCESS = 0; 
  const size_t ERROR_UNHANDLED_EXCEPTION = 2; 
  const size_t HELP_CALLED = 3;
  
  namespace po = boost::program_options; 
 
} /// namespace  

int parseOptionsWithBoost(po::variables_map &vm, po::options_description desc, 
                          int argc, char* argv[]);
TStyle* AtlasStyle();
void SetAtlasStyle();

#endif