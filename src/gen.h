class TRandom3;
class DatabasePDG2;
class Particle;

//#define PLOTS

namespace gen {
// typedef std::vector<Particle*> ParticleList ; // TODO in far future
// data
extern DatabasePDG2 *database;
extern TRandom3 *rnd;

// functions
void load(char *filename, int N);
void initCalc(void);
void doCalculations(int pid = 3122);
void outputPolarization(char *out_file);
void calcInvariantQuantities();
void calcEP1();
}
