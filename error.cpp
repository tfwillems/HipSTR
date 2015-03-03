#include <stdlib.h>
#include <iostream>

#include "error.h"

void printErrorAndDie(std::string message){
  std::cerr << "ERROR: "    << message    << "\n" 
	    << "Exiting..." << std::endl;
  exit(1);
}
