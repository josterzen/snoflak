#include "lattice.h"

#include <math.h>
#include <pthread.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>

using namespace std;

static std::mt19937 generator;

int numOfThreads = 16;

struct threadParameters {
	int start;
	int end;
	Lattice *lattice;
};

Lattice::Lattice(size_t size, std::map<std::string, double> parameters)
	: latticeSize(size), parameters(parameters) {
	resize(size, vector<Cell>(size));
	freezeCenter();
	setDensityOutsideSnowflake();
	updateBoundaryAndComplementOfClosure();
	char *tr = getenv("SLURM_CPUS_PER_TASK");
	numOfThreads = atoi(tr);
}

void Lattice::freezeCenter() {
	size_t middle = ceil(latticeSize / 2);

	(*this)[middle][middle].inSnowflake_before = true;
	(*this)[middle][middle].solidMass_before = 1;
}

void *Lattice::threadSetDensityOutsideSnowflake(void *arg) {
	struct threadParameters pr = *(struct threadParameters *)arg;
	int start = pr.start;
	int end = pr.end;
	size_t stop = end;

	for (size_t i = start; i < stop; i++) {
		for (size_t j = i; j < pr.lattice->latticeSize; j++) {
			(*pr.lattice)[i][j].vaporMass_before =
				pr.lattice->parameters["rho"];
			(*pr.lattice)[j][i].vaporMass_before =
				pr.lattice->parameters["rho"];
		}
	}

	return NULL;
}

void Lattice::setDensityOutsideSnowflake() {
	size_t middle = ceil(latticeSize / 2);

	int currentIndex = 0;
	int end = static_cast<int>(latticeSize);
	int step = end / numOfThreads;

	pthread_t *thread = (pthread_t *)calloc(numOfThreads, sizeof(pthread_t));
	struct threadParameters pr[numOfThreads];

	for (int i = 0; i < numOfThreads - 1; i++) {
		pr[i].start = currentIndex;
		currentIndex += step;
		pr[i].end = currentIndex;
		pr[i].lattice = this;

		pthread_create(&thread[i], NULL,
					   &Lattice::threadSetDensityOutsideSnowflake, &pr[i]);
	}

	pr[numOfThreads - 1].start = currentIndex;
	pr[numOfThreads - 1].end = end;
	pr[numOfThreads - 1].lattice = this;
	pthread_create(&thread[numOfThreads - 1], NULL,
				   &Lattice::threadSetDensityOutsideSnowflake,
				   &pr[numOfThreads - 1]);

	for (int i = 0; i < numOfThreads; i++) {
		pthread_join(thread[i], NULL);
	}

	free(thread);

	(*this)[middle][middle].vaporMass_before = 0;
}

list<pair<long int, long int>> Lattice::getNeighboursIndicesOf(
	long int i, long int j) const {
	int parity = i % 2;

	list<pair<long int, long int>> neighbours{
		{i, j},		{i - 1, j - parity}, {i - 1, j + 1 - parity}, {i, j - 1},
		{i, j + 1}, {i + 1, j - parity}, {i + 1, j + 1 - parity}};

	filterInvalidIndices(neighbours);
	return neighbours;
}

void Lattice::filterInvalidIndices(
	std::list<std::pair<long int, long int>> &indexList) const {
	indexList.remove_if(
		[&](const std::pair<long int, long int> indices) -> bool {
			auto i = indices.first;
			auto j = indices.second;
			return (i < 0) || (j < 0) || (i > latticeSize - 1) ||
				   (j > latticeSize - 1);
		});
}

void Lattice::updateBoundaryAndComplementOfClosure() {
	// printf("Updating boundary and complement of closure");
	// fflush(stdout);

	boundary.clear();
	complementOfClosure.clear();

	std::list<std::pair<long int, long int>> neighbours;
	bool inBoundary;

	for (size_t i = 0; i < latticeSize; i++) {
		for (size_t j = 0; j < latticeSize; j++) {
			inBoundary = false;

			if ((*this)[i][j].inSnowflake_before == true)
				continue;

			neighbours = getNeighboursIndicesOf(i, j);
			// possiblity of using any_of and lambda function
			for (auto &neighbour : neighbours) {
				if ((*this)[neighbour.first][neighbour.second]
						.inSnowflake_before) {
					inBoundary = true;
					break;
				}
			}

			if (inBoundary)
				boundary.push_back({i, j});
			else
				complementOfClosure.push_back({i, j});
		}
	}
	// printf("...done non thread\n");
	// fflush(stdout);
}

void *Lattice::threadDiffuse(void *arg) {
	struct threadParameters pr = *(struct threadParameters *)arg;
	int start = pr.start;
	int end = pr.end;

	list<pair<long int, long int>> neighbours;
	double sum;

	bool inBoundary = false;
	for (int i = start; i < end; i++) {
		sum = 0;

		neighbours = pr.lattice->getNeighboursIndicesOf(
			pr.lattice->complementOfClosure[i].first,
			pr.lattice->complementOfClosure[i].second);

		for (auto &neighbour : neighbours) {
			if ((*pr.lattice)[neighbour.first][neighbour.second]
					.inSnowflake_before)
				sum += (*pr.lattice)[pr.lattice->complementOfClosure[i].first]
									[pr.lattice->complementOfClosure[i].second]
										.vaporMass_before;
			else
				sum += (*pr.lattice)[neighbour.first][neighbour.second]
						   .vaporMass_before;
		}
		sum += (7 - neighbours.size()) *
			   (*pr.lattice)[pr.lattice->complementOfClosure[i].first]
							[pr.lattice->complementOfClosure[i].second]
								.vaporMass_before;

		(*pr.lattice)[pr.lattice->complementOfClosure[i].first]
					 [pr.lattice->complementOfClosure[i].second]
						 .vaporMass_after = sum / 7.0;
	}

	return NULL;
}

void *Lattice::threadDiffuseB(void *arg) {
	struct threadParameters pr = *(struct threadParameters *)arg;
	int start = pr.start;
	int end = pr.end;

	list<pair<long int, long int>> neighbours;
	double sum;

	bool inBoundary = false;
	for (int i = start; i < end; i++) {
		sum = 0;

		neighbours = pr.lattice->getNeighboursIndicesOf(
			pr.lattice->boundary[i].first, pr.lattice->boundary[i].second);

		for (auto &neighbour : neighbours) {
			if ((*pr.lattice)[neighbour.first][neighbour.second]
					.inSnowflake_before)
				sum += (*pr.lattice)[pr.lattice->boundary[i].first]
									[pr.lattice->boundary[i].second]
										.vaporMass_before;
			else
				sum += (*pr.lattice)[neighbour.first][neighbour.second]
						   .vaporMass_before;
		}
		sum += (7 - neighbours.size()) *
			   (*pr.lattice)[pr.lattice->boundary[i].first]
							[pr.lattice->boundary[i].second]
								.vaporMass_before;

		(*pr.lattice)[pr.lattice->boundary[i].first]
					 [pr.lattice->boundary[i].second]
						 .vaporMass_after = sum / 7.0;
	}

	return NULL;
}

void Lattice::diffuse() {
	// printf("Diffusing...\n");
	// fflush(stdout);

	int currentIndex = 0;
	int complementOfClosureSize = static_cast<int>(complementOfClosure.size());
	int boundarySize = static_cast<int>(boundary.size());
	int end = complementOfClosureSize;
	int step = end / numOfThreads;

	pthread_t *thread = (pthread_t *)calloc(numOfThreads, sizeof(pthread_t));
	struct threadParameters pr[numOfThreads];

	for (int i = 0; i < numOfThreads - 1; i++) {
		pr[i].start = currentIndex;
		currentIndex += step;
		pr[i].end = currentIndex;
		pr[i].lattice = this;

		pthread_create(&thread[i], NULL, &Lattice::threadDiffuse, &pr[i]);
	}

	pr[numOfThreads - 1].start = currentIndex;
	pr[numOfThreads - 1].end = end;
	pr[numOfThreads - 1].lattice = this;
	pthread_create(&thread[numOfThreads - 1], NULL, &Lattice::threadDiffuse,
				   &pr[numOfThreads - 1]);

	for (int i = 0; i < numOfThreads; i++) {
		pthread_join(thread[i], NULL);
	}

	currentIndex = 0;
	end = boundarySize;
	step = end / numOfThreads;

	for (int i = 0; i < numOfThreads - 1; i++) {
		pr[i].start = currentIndex;
		currentIndex += step;
		pr[i].end = currentIndex;
		pr[i].lattice = this;

		pthread_create(&thread[i], NULL, &Lattice::threadDiffuseB, &pr[i]);
	}

	pr[numOfThreads - 1].start = currentIndex;
	pr[numOfThreads - 1].end = end;
	pr[numOfThreads - 1].lattice = this;
	pthread_create(&thread[numOfThreads - 1], NULL, &Lattice::threadDiffuseB,
				   &pr[numOfThreads - 1]);

	for (int i = 0; i < numOfThreads; i++) {
		pthread_join(thread[i], NULL);
	}

	free(thread);

	// printf("Diffusion done.\n");
	// fflush(stdout);

	updateValuesOnBoundary();
	updateValuesOnComplementOfClosure();
	updateBoundaryAndComplementOfClosure();
}

void *Lattice::threadFreeze(void *arg) {
	struct threadParameters pr = *(struct threadParameters *)arg;
	int start = pr.start;
	int end = pr.end;

	for (int i = start; i < end; i++) {
		Cell &cell = (*pr.lattice)[pr.lattice->boundary[i].first]
								  [pr.lattice->boundary[i].second];

		cell.liquidMass_after =
			cell.liquidMass_before +
			(1 - pr.lattice->parameters["kappa"]) * cell.vaporMass_before;
		cell.solidMass_after =
			cell.solidMass_before +
			pr.lattice->parameters["kappa"] * cell.vaporMass_before;
		cell.vaporMass_after = 0;
	}

	return NULL;
}

void Lattice::freeze() {
	// printf("Freezing...\n");
	// fflush(stdout);

	int currentIndex = 0;
	int end = boundary.size();
	int step = end / numOfThreads;

	pthread_t *thread = (pthread_t *)calloc(numOfThreads, sizeof(pthread_t));
	struct threadParameters pr[numOfThreads];

	for (int i = 0; i < numOfThreads - 1; i++) {
		pr[i].start = currentIndex;
		currentIndex += step;
		pr[i].end = currentIndex;
		pr[i].lattice = this;

		pthread_create(&thread[i], NULL, &Lattice::threadFreeze, &pr[i]);
	}

	pr[numOfThreads - 1].start = currentIndex;
	pr[numOfThreads - 1].end = end;
	pr[numOfThreads - 1].lattice = this;
	pthread_create(&thread[numOfThreads - 1], NULL, &Lattice::threadFreeze,
				   &pr[numOfThreads - 1]);

	for (int i = 0; i < numOfThreads; i++) {
		pthread_join(thread[i], NULL);
	}

	free(thread);

	updateValuesOnBoundary();
	updateBoundaryAndComplementOfClosure();
}

void *Lattice::threadAttach(void *arg) {
	struct threadParameters pr = *(struct threadParameters *)arg;
	int start = pr.start;
	int end = pr.end;

	list<pair<long int, long int>> neighbours;
	size_t numberOfNeighboursInSnowflake;
	double sumOfVaporMass;

	bool condition1, condition2, condition3;

	for (int i = start; i < end; i++) {
		Cell &cell = (*pr.lattice)[pr.lattice->boundary[i].first]
								  [pr.lattice->boundary[i].second];

		neighbours = pr.lattice->getNeighboursIndicesOf(
			pr.lattice->boundary[i].first, pr.lattice->boundary[i].second);
		neighbours.pop_front();	 // we exclude the cell itself

		numberOfNeighboursInSnowflake = 0;
		sumOfVaporMass = cell.vaporMass_before;
		for (auto &neighbour : neighbours) {
			if ((*pr.lattice)[neighbour.first][neighbour.second]
					.inSnowflake_before)
				numberOfNeighboursInSnowflake++;
			sumOfVaporMass += (*pr.lattice)[neighbour.first][neighbour.second]
								  .vaporMass_before;
		}

		condition1 = (numberOfNeighboursInSnowflake == 1 ||
					  numberOfNeighboursInSnowflake == 2) &&
					 cell.liquidMass_before >= pr.lattice->parameters["beta"];
		condition2 =
			numberOfNeighboursInSnowflake >= 3 &&
			sumOfVaporMass < pr.lattice->parameters["theta"] &&
			cell.liquidMass_before >=
				pr.lattice
					->parameters["alpha"];	// && cell.liquidMass_before >= 1 ;
		condition3 = numberOfNeighboursInSnowflake >= 4;
		if (condition1 || condition2 || condition3) {
			cell.inSnowflake_after = true;
			cell.solidMass_after =
				cell.solidMass_before + cell.liquidMass_before;
			cell.liquidMass_after = 0;
		}
	}
	return NULL;
}

void Lattice::attach() {
	// printf("Attaching...\n");
	// fflush(stdout);

	int currentIndex = 0;
	int end = boundary.size();
	int step = end / numOfThreads;

	pthread_t *thread = (pthread_t *)calloc(numOfThreads, sizeof(pthread_t));
	struct threadParameters pr[numOfThreads];

	for (int i = 0; i < numOfThreads - 1; i++) {
		pr[i].start = currentIndex;
		currentIndex += step;
		pr[i].end = currentIndex;
		pr[i].lattice = this;

		pthread_create(&thread[i], NULL, &Lattice::threadAttach, &pr[i]);
	}

	pr[numOfThreads - 1].start = currentIndex;
	pr[numOfThreads - 1].end = end;
	pr[numOfThreads - 1].lattice = this;
	pthread_create(&thread[numOfThreads - 1], NULL, &Lattice::threadAttach,
				   &pr[numOfThreads - 1]);

	for (int i = 0; i < numOfThreads; i++) {
		pthread_join(thread[i], NULL);
	}

	free(thread);

	updateValuesOnBoundary();
	updateBoundaryAndComplementOfClosure();
}

void *Lattice::threadMelt(void *arg) {
	struct threadParameters pr = *(struct threadParameters *)arg;
	int start = pr.start;
	int end = pr.end;

	for (int i = start; i < end; i++) {
		Cell &cell = (*pr.lattice)[pr.lattice->boundary[i].first]
								  [pr.lattice->boundary[i].second];

		cell.liquidMass_after =
			(1 - pr.lattice->parameters["mu"]) * cell.liquidMass_before;
		cell.solidMass_after =
			(1 - pr.lattice->parameters["gamma"]) * cell.solidMass_before;
		cell.vaporMass_after =
			cell.vaporMass_before +
			pr.lattice->parameters["mu"] * cell.liquidMass_before +
			pr.lattice->parameters["gamma"] * cell.solidMass_before;
	}
	return NULL;
}

void Lattice::melt() {
	// printf("Melting...\n");
	// fflush(stdout);

	int currentIndex = 0;
	int end = boundary.size();
	int step = end / numOfThreads;

	pthread_t *thread = (pthread_t *)calloc(numOfThreads, sizeof(pthread_t));
	struct threadParameters pr[numOfThreads];

	for (int i = 0; i < numOfThreads - 1; i++) {
		pr[i].start = currentIndex;
		currentIndex += step;
		pr[i].end = currentIndex;
		pr[i].lattice = this;

		pthread_create(&thread[i], NULL, &Lattice::threadFreeze, &pr[i]);
	}

	pr[numOfThreads - 1].start = currentIndex;
	pr[numOfThreads - 1].end = end;
	pr[numOfThreads - 1].lattice = this;
	pthread_create(&thread[numOfThreads - 1], NULL, &Lattice::threadFreeze,
				   &pr[numOfThreads - 1]);

	for (int i = 0; i < numOfThreads; i++) {
		pthread_join(thread[i], NULL);
	}

	free(thread);

	updateValuesOnBoundary();
	updateBoundaryAndComplementOfClosure();
}

void *Lattice::threadAddNoise(void *arg) {
	struct threadParameters pr = *(struct threadParameters *)arg;
	int start = pr.start;
	int end = pr.end;

	std::bernoulli_distribution bernoulliDistribution(0.5);
	double perturbation;
	int sum = 0;

	bool inBoundary = false;
	for (int i = start; i < end; i++) {
		sum = 0;

		Cell &cell = (*pr.lattice)[pr.lattice->complementOfClosure[i].first]
								  [pr.lattice->complementOfClosure[i].second];

		perturbation = (2 * bernoulliDistribution(generator) - 1) *
					   pr.lattice->parameters["sigma"];
		cell.vaporMass_after = (1 - perturbation) * cell.vaporMass_before;
	}

	return NULL;
}

void *Lattice::threadAddNoiseB(void *arg) {
	struct threadParameters pr = *(struct threadParameters *)arg;
	int start = pr.start;
	int end = pr.end;

	std::bernoulli_distribution bernoulliDistribution(0.5);
	double perturbation;
	int sum = 0;

	bool inBoundary = false;
	for (int i = start; i < end; i++) {
		sum = 0;

		Cell &cell = (*pr.lattice)[pr.lattice->boundary[i].first]
								  [pr.lattice->boundary[i].second];

		perturbation = (2 * bernoulliDistribution(generator) - 1) *
					   pr.lattice->parameters["sigma"];
		cell.vaporMass_after = (1 - perturbation) * cell.vaporMass_before;
	}

	return NULL;
}

void Lattice::addNoise() {
	// printf("Adding noise...\n");
	// fflush(stdout);

	int currentIndex = 0;
	int complementOfClosureSize = static_cast<int>(complementOfClosure.size());
	int boundarySize = static_cast<int>(boundary.size());
	int end = complementOfClosureSize;
	int step = end / numOfThreads;

	pthread_t *thread = (pthread_t *)calloc(numOfThreads, sizeof(pthread_t));
	struct threadParameters pr[numOfThreads];

	for (int i = 0; i < numOfThreads - 1; i++) {
		pr[i].start = currentIndex;
		currentIndex += step;
		pr[i].end = currentIndex;
		pr[i].lattice = this;

		pthread_create(&thread[i], NULL, &Lattice::threadAddNoise, &pr[i]);
	}

	pr[numOfThreads - 1].start = currentIndex;
	pr[numOfThreads - 1].end = end;
	pr[numOfThreads - 1].lattice = this;
	pthread_create(&thread[numOfThreads - 1], NULL, &Lattice::threadAddNoise,
				   &pr[numOfThreads - 1]);

	for (int i = 0; i < numOfThreads; i++) {
		pthread_join(thread[i], NULL);
	}

	currentIndex = 0;
	end = boundarySize;
	step = end / numOfThreads;

	for (int i = 0; i < numOfThreads - 1; i++) {
		pr[i].start = currentIndex;
		currentIndex += step;
		pr[i].end = currentIndex;
		pr[i].lattice = this;

		pthread_create(&thread[i], NULL, &Lattice::threadAddNoiseB, &pr[i]);
	}

	pr[numOfThreads - 1].start = currentIndex;
	pr[numOfThreads - 1].end = end;
	pr[numOfThreads - 1].lattice = this;
	pthread_create(&thread[numOfThreads - 1], NULL, &Lattice::threadAddNoiseB,
				   &pr[numOfThreads - 1]);

	for (int i = 0; i < numOfThreads; i++) {
		pthread_join(thread[i], NULL);
	}

	free(thread);

	updateValuesOnBoundary();
	updateValuesOnComplementOfClosure();
	updateBoundaryAndComplementOfClosure();
}

void *Lattice::threadUpdateValuesOnBoundary(void *arg) {
	struct threadParameters pr = *(struct threadParameters *)arg;
	int start = pr.start;
	int end = pr.end;

	for (int i = start; i < end; i++) {
		Cell &cell = (*pr.lattice)[pr.lattice->boundary[i].first]
								  [pr.lattice->boundary[i].second];

		cell.liquidMass_before = cell.liquidMass_after;
		cell.solidMass_before = cell.solidMass_after;
		cell.vaporMass_before = cell.vaporMass_after;
		cell.inSnowflake_before = cell.inSnowflake_after;
	}
	return NULL;
}

void Lattice::updateValuesOnBoundary() {
	// printf("Updating values on boundary...\n");
	// fflush(stdout);

	int currentIndex = 0;
	int end = boundary.size();
	int step = end / numOfThreads;

	pthread_t *thread = (pthread_t *)calloc(numOfThreads, sizeof(pthread_t));
	struct threadParameters pr[numOfThreads];

	for (int i = 0; i < numOfThreads - 1; i++) {
		pr[i].start = currentIndex;
		currentIndex += step;
		pr[i].end = currentIndex;
		pr[i].lattice = this;

		pthread_create(&thread[i], NULL, &Lattice::threadUpdateValuesOnBoundary,
					   &pr[i]);
	}

	pr[numOfThreads - 1].start = currentIndex;
	pr[numOfThreads - 1].end = end;
	pr[numOfThreads - 1].lattice = this;
	pthread_create(&thread[numOfThreads - 1], NULL,
				   &Lattice::threadUpdateValuesOnBoundary,
				   &pr[numOfThreads - 1]);

	for (int i = 0; i < numOfThreads; i++) {
		pthread_join(thread[i], NULL);
	}

	free(thread);

	// printf("Done updating values on boundary.\n");
	// fflush(stdout);
}

void *Lattice::threadUpdateValuesOnComplementOfClosure(void *arg) {
	// printf("Thread start\n");
	// fflush(stdout);

	struct threadParameters pr = *(struct threadParameters *)arg;
	int start = pr.start;
	int end = pr.end;

	for (int i = start; i < end; i++) {
		Cell &cell = (*pr.lattice)[pr.lattice->complementOfClosure[i].first]
								  [pr.lattice->complementOfClosure[i].second];

		cell.liquidMass_before = cell.liquidMass_after;
		cell.solidMass_before = cell.solidMass_after;
		cell.vaporMass_before = cell.vaporMass_after;
		cell.inSnowflake_before = cell.inSnowflake_after;
	}

	// printf("Thread end\n");
	// fflush(stdout);
	return NULL;
}

void Lattice::updateValuesOnComplementOfClosure() {
	// printf("Updating values on complement of closure...\n");
	// fflush(stdout);

	int currentIndex = 0;
	int end = complementOfClosure.size();
	int step = end / numOfThreads;

	pthread_t *thread = (pthread_t *)calloc(numOfThreads, sizeof(pthread_t));
	struct threadParameters pr[numOfThreads];

	for (int i = 0; i < numOfThreads - 1; i++) {
		pr[i].start = currentIndex;
		currentIndex += step;
		pr[i].end = currentIndex;
		pr[i].lattice = this;

		pthread_create(&thread[i], NULL,
					   &Lattice::threadUpdateValuesOnComplementOfClosure,
					   &pr[i]);
	}

	pr[numOfThreads - 1].start = currentIndex;
	pr[numOfThreads - 1].end = end;
	pr[numOfThreads - 1].lattice = this;
	pthread_create(&thread[numOfThreads - 1], NULL,
				   &Lattice::threadUpdateValuesOnComplementOfClosure,
				   &pr[numOfThreads - 1]);

	for (int i = 0; i < numOfThreads; i++) {
		pthread_join(thread[i], NULL);
	}

	free(thread);

	// printf("Done updating values on complement of closure.\n");
	// fflush(stdout);
}

void Lattice::saveTo(std::string fileName) const {
	ofstream inSnowflake_fileStream(fileName + "_inSnowflake",
									ios::out | ios::binary);
	ofstream vaporMass_fileStream(fileName + "_vapor", ios::out | ios::binary);
	ofstream liquidMass_fileStream(fileName + "_liquid",
								   ios::out | ios::binary);
	ofstream solidMass_fileStream(fileName + "_solid", ios::out | ios::binary);

	if (!inSnowflake_fileStream.is_open() || !vaporMass_fileStream.is_open() ||
		!liquidMass_fileStream.is_open() || !solidMass_fileStream.is_open())
		throw runtime_error(
			"Could not open the file to export the lattice. Make sure the "
			"output folder was created.");

	for (size_t i = 0; i < latticeSize; i++) {
		for (size_t j = 0; j < latticeSize; j++) {
			const Cell &cell = (*this)[i][j];

			inSnowflake_fileStream.write((char *)&(cell.inSnowflake_before),
										 sizeof(cell.inSnowflake_before));
			liquidMass_fileStream.write((char *)&(cell.liquidMass_before),
										sizeof(cell.liquidMass_before));
			vaporMass_fileStream.write((char *)&(cell.vaporMass_before),
									   sizeof(cell.vaporMass_before));
			solidMass_fileStream.write((char *)&(cell.solidMass_before),
									   sizeof(cell.solidMass_before));
		}
	}
}
