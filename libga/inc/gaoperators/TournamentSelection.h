#ifndef __TOURNAMENTSELECTION__
#define __TOURNAMENTSELECTION__

#include "Selection.h"

class TournamentSelection : public Selection<TournamentSelection> {

public:
  template <typename F, typename T>
  static void
  SelectionGA(Population<T> &oldpop,
              Population<T> &newpop) throw(ExceptionMessenger) {
    std::random_device rd;
    std::mt19937 gen(rd());
    const Int_t PopSizeCheck = oldpop.GetPopulationSize();
    if ((newpop.GetPopulationSize()) != PopSizeCheck)
      throw ExceptionMessenger("OMG! New population has wrong size");
    std::vector<Int_t> Population1(PopSizeCheck), Population2(PopSizeCheck);
    for (Int_t i = 0; i < PopSizeCheck; ++i) {
      Population1[i] = Population2[i] = i;
    }
    for (Int_t i = 0; i < PopSizeCheck; ++i) {
      std::uniform_int_distribution<> dist(i, PopSizeCheck - 1);
      // std::swap(Population1[rand.Integer(PopSizeCheck - 1)], Population1[i]);
      // std::swap(Population2[rand.Integer(PopSizeCheck - 1)], Population2[i]);
      // Int_t rand1 = rand() % (PopSizeCheck - 1) + i;
      // std::swap(Population1[rand1], Population1[i]);
      // Int_t rand2 = rand() % (PopSizeCheck - 1) + i;
      // std::swap(Population2[rand2], Population2[i]);
      std::swap(Population1[dist(gen)], Population1[i]);
      std::swap(Population2[dist(gen)], Population2[i]);
    }
    for (Int_t i = 0; i < PopSizeCheck; i += 4) {
      Genes<T> &Combination11 = Tournament(
          oldpop.GetGenes(Population1[i]), oldpop.GetGenes(Population1[i + 1]));
      Genes<T> &Combination12 =
          Tournament(oldpop.GetGenes(Population1[i + 2]),
                     oldpop.GetGenes(Population1[i + 3]));
      // std::cout << "-================== Combination11 & Combination12
      // ============================-" << std::endl;
      // Combination11.Genes<T>::printGenes(Combination11);
      // std::cout << "-==============================================-" <<
      // std::endl;
      // Combination12.Genes<T>::printGenes(Combination12);
      Crossover(Combination11, Combination12, newpop.GetGenes(i),
                newpop.GetGenes(i + 1));
      Genes<T> &Combination21 = Tournament(
          oldpop.GetGenes(Population2[i]), oldpop.GetGenes(Population2[i + 1]));
      Genes<T> &Combination22 =
          Tournament(oldpop.GetGenes(Population2[i + 2]),
                     oldpop.GetGenes(Population2[i + 3]));
      // std::cout << "-================== Combination21 & Combination22
      // ============================-" << std::endl;
      // Combination21.Genes<T>::printGenes(Combination21);
      // std::cout << "-==============================================-" <<
      // std::endl;
      // Combination22.Genes<T>::printGenes(Combination22);
      // std::cout << "-==============================================-" <<
      // std::endl;
      Crossover(Combination21, Combination22, newpop.GetGenes(i + 2),
                newpop.GetGenes(i + 3));
    }
  }
};

#endif
