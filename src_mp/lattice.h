#ifndef LATTICE_H
#define LATTICE_H

#include <list>
#include <map>
#include <string>
#include <vector>

struct Cell {
	bool inSnowflake_before = false;
	double vaporMass_before = 0;
	double solidMass_before = 0;
	double liquidMass_before = 0;

	bool inSnowflake_after = false;
	double vaporMass_after = 0;
	double solidMass_after = 0;
	double liquidMass_after = 0;
};

class Lattice : public std::vector<std::vector<Cell>> {
	size_t latticeSize;
	std::map<std::string, double> parameters;

	std::vector<std::pair<size_t, size_t>> boundary;
	std::vector<std::pair<size_t, size_t>> complementOfClosure;

   public:
	Lattice(size_t size, std::map<std::string, double> parameters);

	void diffuse();
	void freeze();
	void attach();
	void melt();
	void addNoise();

	void saveTo(std::string fileName) const;

   private:
	void freezeCenter();
	void setDensityOutsideSnowflake();

	void updateBoundaryAndComplementOfClosure();
	void updateValuesOnBoundary();
	void updateValuesOnComplementOfClosure();

	std::list<std::pair<long int, long int>> getNeighboursIndicesOf(
		long int i, long int j) const;
	void filterInvalidIndices(
		std::list<std::pair<long int, long int>> &indexList) const;
};

#endif
