#ifndef ROOTS_H
#define ROOTS_H

#include <string>
#include <utility>
#include <random>
#include <cstddef>

class roots
/***
  Ineficient but fun way to find roots to continuous functions.
***/
{
private:
  std::string* genes;
  std::size_t numGenes;
  int decimals, crossPoints;
  double lower, upper, mutRate;
  std::default_random_engine rnd;
  double (*f)(double); //the function which we want to find roots of
  double fitness(std::string);
  std::string tostring(double);
  void mutateGene(int);
  void mutate();
  std::string crossnums(std::string,std::string, std::vector<std::size_t>);
  std::vector<std::size_t> makePoints(std::size_t);
  std::string crossover(std::string,std::string);
  void repopulate(std::string,std::string);

public:
  roots(std::size_t, double, double, double(*) (double));
  roots(const roots&);
  roots& operator=(roots&);
  ~roots();
  void setMutation(double);
  void setF(double(*) (double));
  void setNumGenes(std::size_t);
  void initialize();
  std::pair<std::size_t,std::size_t> top_two();
  void generation(); //equivalent to run(1)
  void run(int);
  std::string fittest();

  friend std::ostream& operator << (std::ostream&, const roots);
};



















#endif
