class TRandom3;
class DatabasePDG2;
class Particle;

namespace gen {
// typedef std::vector<Particle*> ParticleList ; // TODO in far future
// data
extern DatabasePDG2 *database;
extern TRandom3 *rnd;

// functions
void load(char *filename, int N);
void initCalc(void);
void doCalculations(void);
void calcBasicQuantities(void);
}
