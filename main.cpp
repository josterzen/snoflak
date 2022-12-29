#include <math.h>

#include <fstream>
#include <iostream>

#include "lattice.h"

using namespace std;

/*static char* getCharFromString(const string& str) {
	char* charPtr = new char[str.length() + 1];
	strcpy(charPtr, str.c_str());
	return charPtr;
}

static rapidjson::Document getJsonDocumentFromFile(const string& fileName) {
	ifstream fileStream(fileName);
	if (!fileStream.is_open())
		throw runtime_error("Could not open config file");
	std::string tmp((istreambuf_iterator<char>(fileStream)),
					istreambuf_iterator<char>());
	const char* fileContent = getCharFromString(tmp);

	rapidjson::Document document;
	document.Parse<rapidjson::kParseCommentsFlag>(fileContent);
	delete[] fileContent;
	return document;
}

static std::map<std::string, double> getMapFromJsonDocument(
	rapidjson::Document& document) {
	std::map<std::string, double> parametersMap;

	std::list<std::string> parameters{"rho",   "beta", "alpha", "theta",
									  "kappa", "mu",   "gamma", "sigma"};

	for (auto& parameter : parameters) {
		const char* parameterName = getCharFromString(parameter);
		parametersMap.emplace(parameter, document[parameterName].GetDouble());
		delete[] parameterName;
	}
	return parametersMap;
}*/

int main(int argc, char* argv[]) {
	if (argc == 1) {
		std::cout << "You must enter the json file containing the parameters"
				  << std::endl;
		return 1;
	}

	bool skipFrames = true;

	// auto document = getJsonDocumentFromFile(argv[1]);
	std::map<std::string, double> parameters = {
		{"size", 600}, {"iterations", 100}, {"frameskips", 0}, {"rho", 0.58},
		{"beta", 2},   {"alpha", 0.8},		{"theta", 0.011},  {"kappa", 0.1},
		{"mu", 0.01},  {"gamma", 0.00005},	{"sigma", 0.00005}};

	size_t iterations = 1000;  // document["iterations"].GetInt();
	size_t size = 600;

	Lattice lattice(size, parameters);

	for (size_t i = 0; i < iterations; i++) {
		std::cout << "iterations " + std::to_string(i) + "/" +
						 std::to_string(iterations)
				  << "\r";
		lattice.diffuse();
		lattice.freeze();
		lattice.attach();
		lattice.melt();
		lattice.addNoise();

		if (skipFrames) {
			lattice.saveTo("frames/snowflake" + std::string("_") +
						   std::to_string(i));
		}
	}
	lattice.saveTo("frames/snowflake");

	return 0;
}
