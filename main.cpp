#include "roots.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <typeinfo>
#include <utility>
#include <unordered_set>
#include <random>
#include <vector>
#include <cstddef>

double foo(double x)
{
  return x*x*x + x*x - 6*x;
}



int main()
{

  roots R(100,1,3,foo);
  R.initialize();
  R.setMutation(0.2);
  //std::cout<<R<<std::endl;
  R.run(500);
  std::cout<<R.fittest()<<std::endl;

  return 0;
}
