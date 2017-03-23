#include <iostream>
#include <math.h>

#include "../src/fastonebigheader.h"

int main(){
  for (double v = -20; v < 0; v += 0.1){
    std::cerr << v << "\t" << exp(v) << "\t" << fastexp(v) << std::endl;
  }

  for (double v = -20; v < 0; v += 0.1){
    double x = exp(v);
    std::cerr << v << "\t" << log(x) << "\t" << fastlog(x) << std::endl;
  }
}
