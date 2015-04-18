//---------------------------------------------------------------------------
//
// Name: Philippe Demontigny
//
// Student Number: 20557658
// User-id: pdemonti
// Assignment: 4
//
//---------------------------------------------------------------------------


#include <iostream>
#include "scene_lua.hpp"

bool bounding = false;
bool fast = true;

int main(int argc, char** argv)
{
  std::string filename = "scene.lua";
  if (argc >= 3) {
  	std::string cond = argv[1];
  	if (cond == "slow") {
  		fast = false;
  	}
  	else if (cond == "fast") {
  		fast = true;
  	}
  	else if (cond == "bounding") {
  		bounding = true;
  	}
  	else {
  		std::cerr << "ERROR: Command " << cond << " not recognized." << std::endl;
  	}
    filename = argv[2];
  }
  else if (argc >= 2) {
  	filename = argv[1];
  }

  if (!run_lua(filename)) {
    std::cerr << "Could not open " << filename << std::endl;
    return 1;
  }
}

