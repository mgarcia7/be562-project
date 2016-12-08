#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>


int gap_penalty = 4;
int S[4][4] = { {3, -1, -2, -2},
				{-1, 3, -2, -2},
			    {-2, -2, 3, -1},
				{-2, -2, -1, 3} };
std::map<char, int> base_idx = {{'A', 0}, {'G', 1}, {'C', 2}, {'T', 3}};
int PTR_NONE = 0;
int PTR_GAP1 = 1;
int PTR_GAP2 = 2;
int PTR_BASE = 3;


std::string read_fasta(std::string file_name) {
	std::ifstream handle(file_name);
	std::string l;
	std::string genome;
	int count = 1;

	while (std::getline(handle, l)) {
		if (count != 1) {
			//skip first line starting with '>'
			continue;
		} else {
			genome.append(l);
		}
	}

	return genome;
}


void Smith_Waterman(std::string& s1, std::string& s2,
					std::vector<std::vector<int> >* F,
					std::vector<std::vector<int> >* TB) {

	for (std::size_t i = 1; i < F->size(); i++) {
		for (std::size_t j = 1; j < (*F)[i].size(); j++) {

			std::vector<int> scores(3, 0);
			scores[0] = (*F)[i-1][j-1] + S[base_idx[s1[i-1]]][base_idx[s2[j-1]]];
			scores[1] = (*F)[i][j-1] - gap_penalty;
			scores[2] = (*F)[i-1][j] - gap_penalty;

			int max = 0;
			int max_idx = 0;
			for (std::size_t a = 0; a < 3; a++) {
				if (scores[a] > max) {
					max = scores[a];
					max_idx = a;
				}
			}

			if (max_idx == 0) {
				(*TB)[i][j] = PTR_BASE;
			} else if (max_idx == 1) {
				(*TB)[i][j] = PTR_GAP1;
			} else {
				(*TB)[i][j] = PTR_GAP2;
			}

			if (max < 0) {
				(*F)[i][j] = 0;
			} else {
				(*F)[i][j] = max;
			}

		}
	}

}


void find_max(std::vector<std::vector<int> >* M, std::size_t& i, std::size_t& j) {
	int t = 0;
	std::size_t i_max = 0;
	std::size_t j_max = 0;
	for (std::size_t i = 0; i < (*M).size(); i++) {
		for (std::size_t j = 0; j < (*M)[i].size(); j++) {
			if ((*M)[i][j] > t) {
				t = (*M)[i][j];
				i_max = i;
				j_max = i;
			}
		}
	}

	i = i_max;
	j = j_max;

}


void traceback(std::string& s1, std::string& s2,
			   std::vector<std::vector<int> >* F,
			   std::vector<std::vector<int> >* TB,
			   std::vector<std::string>& conserved) {

	int max_iter = 5;
	int iterations = 0;

	std::string c1("");
	std::string c2("");

	std::size_t i = 0;
	std::size_t j = 0;
	find_max(F, i, j);

	while (iterations < max_iter) {
		iterations += 1;
		while ((*F)[i][j] != 0) {
			if ((*TB)[i][j] == PTR_BASE) {
				(*F)[i][j] = 0;
				c1 = s1[i-1] + c1;
				c2 = s2[j-1] + c2;
				i = i - 1;
				j = j - 1;
			} else if ((*TB)[i][j] == PTR_GAP1) {
				(*F)[i][j] = 0;
				c1 = '-' + c1;
				c2 = s2[j-1] + c2;
				j = j - 1;
			} else if ((*TB)[i][j] == PTR_GAP2) {
				(*F)[i][j] = 0;
				c1 = s1[i-1] + c1;
				c2 = '-' + c2;
				i = i - 1;
			}
		}

		conserved.push_back(c1);
		conserved.push_back(c2);
		c1.clear();
		c2.clear();
		i = 0;
		j = 0;
		find_max(F, i, j);
	}

}


int main(int argc, const char* argv[]) {

	//convert input args to strings
	std::vector<std::string> file_names;
	for (int i = 1; i < argc; i++) {
		file_names.push_back(argv[i]);
	}

	//read genomes into vector of strings
	std::vector<std::string> genomes;
	for (std::string f : file_names) {
		genomes.push_back(read_fasta(f));
	}

	auto F = new std::vector<std::vector<int> >(genomes[0].size()+1, std::vector<int>(genomes[1].size()+1, 0));
	auto TB = new std::vector<std::vector<int> >(genomes[0].size()+1, std::vector<int>(genomes[1].size()+1, PTR_NONE));

	Smith_Waterman(genomes[0], genomes[1], F, TB);

	std::vector<std::string> conserved;
	traceback(genomes[1], genomes[1], F, TB, conserved);

	for (std::string s : conserved) {
		std::cout << s << std::endl << std::endl;
	}

	delete F;
	delete TB;

	return 0;

}