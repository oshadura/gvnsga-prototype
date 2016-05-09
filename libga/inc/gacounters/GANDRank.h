#ifndef __GANDRANK__
#define __GANDRANK__

class GANDRank : public GACounter<GANDRank, double> {
public:
  template <typename F>
      static std::unordered_map < std::vector<std::shared_ptr<Genes<F> > >
      CalculateRank(Population<F> pop, int fBestIndividuals = -1) {
    if (fBestIndividuals == -1)
      fBestIndividual = pop.size();
    std::unordered_map<std::vector<std::shared_ptr<Genes<F> > >, int> fMap;
    int fRank = 0;
    int fNumberOfIndividuals = 0;
    while (!pop.empty() && (fNumberOfIndividuals < fBestIndividuals)) {
      auto fFront = UpdatePF::GetParetoFront(pop);
      pop.remove(fFront);
      for (std::vector<std::shared_pointers<Genes<F>>> entry : fFront) {
        fMap[entry] = fRank;
        ++fNumberOfIndividuals;
      }
      ++fRank;
    }
    for (std::vector<std::shared_ptr<Genes<F> > > entry : pop)
      fMap[entry] = std::numeric_limits<int>::max();
    return fMap;
  }
};

#endif