#include "roots.h"
#include <string>
#include <cstdlib>
#include <stdlib.h>
#include <iostream>
#include <utility>
#include <sstream>
#include <vector>
#include <algorithm>
#include <random>
#include <cstddef>
#include <unordered_set>

roots::roots(std::size_t n, double l, double u, double(*_f) (double))
{
  numGenes = n;
  genes = new std::string[numGenes];
  lower = l;
  upper = u;
  f=_f;
  decimals = 40;
  crossPoints = 20;
}

roots::roots(const roots& oldroots)
{
  numGenes = oldroots.numGenes;
  genes = new std::string[numGenes];
  for (std::size_t i=0;i<numGenes;i++)
  {
    genes[i] = oldroots.genes[i];
  }
  lower=oldroots.lower;
  upper=oldroots.upper;
  f=oldroots.f;
  decimals = oldroots.decimals;
  crossPoints = oldroots.crossPoints;
}

roots& roots::operator = (roots& oldroots)
{
  numGenes = oldroots.numGenes;
  for (std::size_t i=0;i<numGenes;i++)
  {
    genes[i] = oldroots.genes[i];
  }
  lower=oldroots.lower;
  upper=oldroots.upper;
  f=oldroots.f;
  decimals = oldroots.decimals;
}

roots::~roots()
{
  if (genes)
  {
    delete[] genes;
  }
}

void roots::setMutation(double newMut)
{
  mutRate = newMut;
}

void roots::setF(double(*_f) (double))
{
  f=_f;
}

void roots::setNumGenes(std::size_t n)
{
  numGenes = n;
}

double roots::fitness(std::string gene)
{
  try
  {
    double intgene = std::stod(gene);
    return std::abs(f(intgene));
  }
  catch (std::invalid_argument)
  {
    std::cout<<"bad gene: "<<gene<<std::endl;
  }

}

std::string roots::tostring(double x)
{
    std::ostringstream out;
    out.precision(decimals);
    out<<std::fixed<<x;
    std::string result=out.str();
    return result;
}

void roots::initialize()
{
  std::uniform_real_distribution<double> distribution(lower,upper);
  for (std::size_t i=0;i<numGenes;i++)
  {
    genes[i] = tostring(distribution(rnd));
  }
}

void roots::mutateGene(int index)
/***
***/
{
  std::string gene = genes[index];
  std::vector<char> geneVec(gene.begin(),gene.end());
  std::size_t dot = gene.find_first_of('.');
  std::vector<char>::iterator zero = geneVec.begin();
  geneVec.erase(zero+dot);

  bool sign = (gene.at(0)=='-');
  if (sign)
  {
    geneVec.erase(zero);
  }
  for (int i=0;i<geneVec.size();i++)
  {
    std::uniform_real_distribution<double> shouldMut(0,1);
    if (shouldMut(rnd) < mutRate)
    {
      std::uniform_int_distribution<int> digit(0,9);
      geneVec.at(i) = 48+digit(rnd);
    }

  }
  if (sign)
  {
    geneVec.insert(zero,'-');
  }
  geneVec.insert(zero+dot,'.');
  genes[index]=std::string(zero,geneVec.end());
}

void roots::mutate()
{
  for (std::size_t i=0;i<numGenes;i++)
  {
    mutateGene(i);
  }
}

std::pair<std::size_t,std::size_t> roots::top_two()
/***
Returns the fittest two genes in the gene pool.
***/
{

  std::size_t first = 0;
  double ffitness = fitness(genes[0]);
  std::size_t second = 1;
  double sfitness = fitness(genes[1]);
  for (std::size_t i=1;i<numGenes;i++)
  {
    double ifitness = fitness(genes[i]);
    if (ifitness < ffitness)
    {
      second = first;
      first=i;
      sfitness = ffitness;
      ffitness = ifitness;
    }
    else if (ifitness < sfitness)
    {
      second = i;
      sfitness = ifitness;
    }
  }
  std::pair<std::size_t,std::size_t> best = std::pair<std::size_t,std::size_t>(first, second);
  return best;
}

std::string roots::crossnums(std::string mom, std::string dad, std::vector<std::size_t> points)
{

  int limit = mom.size();
  points.push_back(limit);

  bool toggle=true;
  std::vector<char> preresult;
  std::size_t j=0;
  for (std::size_t i=0;i<limit;i++)
  {
    if (toggle || i>=dad.size())
    {
      preresult.push_back(mom.at(i));
    }
    else
    {
      preresult.push_back(dad.at(i));
    }
    if (i==points.at(j))
    {
      toggle = !toggle;
      j++;
    }
  }
  std::string result = std::string(preresult.begin(), preresult.end());
  return result;
}

std::vector<std::size_t> roots::makePoints(std::size_t len)
{
  std::vector<std::size_t> result;
  std::size_t prev = len+1;
  for (std::size_t i=0;i<crossPoints;i++)
  {
    std::uniform_int_distribution<std::size_t> dist(0,prev-1);
    result.push_back(dist(rnd));
    if (prev==0)
    {
      break;
    }
    prev = result.at(i);
  }
  return result;
}

std::size_t _find(std::vector<std::size_t> arr, std::size_t key, std::size_t l, std::size_t h)
{
  /***
  Finds the location of a key in an ordered vector arr between l and h. If key is less than l, it returns l. If it is higher than h, it returns h+1
  If key is between two elements x and y of arr, it returns the location of y.
  ***/

  std::size_t x = (l+h)/2;
  if (arr.at(x) == key)
  {
    return x;
  }
  if (key <= arr.at(l))
  {
    return l;
  }
  if (key >= arr.at(h-1))
  {
    return h;
  }
  if (h==l+1)
  {
    return h;
  }
  if (key>arr.at(x))
  {
    return _find(arr, key, x, h);
  }
  //else
  return _find(arr, key, l, x);
}

std::size_t find(std::vector<std::size_t> arr, std::size_t key)
{
  return _find(arr, key, 0, arr.size());
}

std::string roots::crossover(std::string mom, std::string dad)
{
  //create sign for the beginning of the gene:
  std::string sign("");
  if (mom.at(0)=='-')
  {
    mom = mom.substr(1);
    if (dad.at(0)=='-')
    {
      sign = "-";
      dad = dad.substr(1);
    }
    else
    {
      std::uniform_int_distribution<int> choose(0,1);
      if (choose(rnd)==0)
      {
        sign="-";
      }
    }
  }
  else if (dad.at(0)=='-')
  {
    dad = dad.substr(1);
    std::uniform_int_distribution<int> choose(0,1);
    if (choose(rnd)==0)
    {
      sign="-";
    }
  }
  //___________

  std::size_t momDot = mom.find_first_of('.');
  std::size_t dadDot = dad.find_first_of('.');
  std::vector<std::size_t> points = makePoints(mom.size()-1);
  std::size_t pointDot = find(points, momDot);
  std::vector<std::size_t> intpoints(points.begin(), points.begin()+pointDot);

  std::string integral = crossnums(mom.substr(0,momDot),dad.substr(0,dadDot),intpoints);
  std::vector<std::size_t> decpoints(points.begin()+pointDot, points.end());
  std::string dec = crossnums(mom.substr(momDot+1), dad.substr(dadDot+1),decpoints);
  return sign+integral+"."+dec;
}



void roots::repopulate(std::string mom, std::string dad)
{
  for (int i=0;i<numGenes;i++)
  {
    if (i%2==0)
    {
      genes[i] = crossover(mom,dad);
    }
    else
    {
      genes[i] = crossover(dad,mom);
    }
  }
}

void roots::generation()
{
  std::pair<std::size_t,std::size_t> best = top_two();
  repopulate(genes[best.first],genes[best.second]);
  mutate();
}

void roots::run(int n)
{
  for (int i=0;i<n;i++)
  {
    generation();
  }
}

std::string roots::fittest()
{
  std::size_t first = 0;
  double ffitness = fitness(genes[0]);
  for (std::size_t i=1;i<numGenes;i++)
  {
    double ifitness = fitness(genes[i]);
    if (ifitness < ffitness)
    {
      first=i;
      ffitness = ifitness;
    }
  }
  return genes[first];
}



std::ostream& operator << (std::ostream& os, const roots R)
{
  for (int i=0;i<R.numGenes;i++)
  {
    std::string gene = R.genes[i];
    os<<gene<<std::endl;
  }
  return os;
}
