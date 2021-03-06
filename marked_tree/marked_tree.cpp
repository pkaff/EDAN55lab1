#include "stdafx.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>

using namespace std;

int N;
#define MAX_LAYERS 20
#define N_ROUNDS 10.0

class Node {
public:

	int id;
	int left;
	int right;
	int parent;
	int layer;

	bool marked = false;

	Node() {}

	Node(int id) : id(id) {
		left = 2 * id;
		right = 2 * id + 1;
		parent = id / 2;
		layer = floor(log2(id)) + 1;
	}
};

bool mark(int i, vector<Node>& graph, bool checkFinish) {
	graph[i].marked = true;

	int finalLayer = graph.back().layer;

	//Rule 1
	if (graph[i].id != 1) {
		if (graph[i].id % 2 == 0) { //left
			if (graph[i + 1].marked) {
				if (!graph[graph[i].parent].marked)
					mark(graph[i].parent, graph, false);
			}
		}
		else if (graph[i].id % 2 == 1) { //right
			if (graph[i - 1].marked) {
				if (!graph[graph[i].parent].marked)
					mark(graph[i].parent, graph, false);
			}
		}
	}

	//Rule 2
	if (graph[i].id != 1) {
		if (graph[graph[i].parent].marked) {
			if (graph[i].id % 2 == 0) { //left
				if (!graph[i + 1].marked)
					mark(i + 1, graph, false);
			}
			else if (graph[i].id % 2 == 1) { //right
				if (!graph[i - 1].marked)
					mark(i - 1, graph, false);
			}
		}
	}
	if (graph[i].layer != finalLayer) {
		if (graph[graph[i].left].marked) {
			if (!graph[graph[i].right].marked)
				mark(graph[i].right, graph, false);
		}
		else if (graph[graph[i].right].marked) {
			if (!graph[graph[i].left].marked)
				mark(graph[i].left, graph, false);
		}
	}

	if (checkFinish) {
		//check for termination
		size_t j = N;
		while (graph[j].layer == finalLayer) {
			if (!graph[j--].marked) {
				return false;
			}
		}
		return true;
	}
	else {
		return false;
	}
}

void resetGraph(vector<Node>& graph) {
	for (int i = 1; i <= N; ++i) {
		graph[i].marked = false;
	}
}

void assertCompletion(int t, const vector<Node>& graph) {
	for (int i = 1; i < graph.size(); ++i) {
		if (!graph[i].marked) {
			cout << "Test " << t << " not terminated correctly!" << endl;
			return;
		}
	}
}

int main() {
	ofstream ofs("C:\\Users\\biz\\Desktop\\output.txt");
	vector<vector<double>> iterations1(N_ROUNDS, vector<double>(MAX_LAYERS - 1, 0.0));
	vector<vector<double>> iterations2(N_ROUNDS, vector<double>(MAX_LAYERS - 1, 0.0));
	vector<vector<double>> iterations3(N_ROUNDS, vector<double>(MAX_LAYERS - 1, 0.0));
	for (int r = 0; r < N_ROUNDS; ++r) {
		for (int l = 2; l <= MAX_LAYERS; ++l) {
			N = pow(2, l) - 1;
			int itr1 = 0;
			int itr2 = 0;
			int itr3 = 0;
			vector<Node> graph(N + 1);
			vector<int> v2(N);
			vector<int> v3(N);
			for (int i = 1; i <= N; ++i) {
				Node n(i);
				graph[i] = n;
				v2[i - 1] = i;
				v3[i - 1] = i;
			}

			random_device rd;
			mt19937 gen(rd());
			uniform_int_distribution<unsigned long long> dis(1, N);
			while (true) {
				unsigned long long randomNbr = dis(gen);
				++itr1;
				if (graph[randomNbr].marked) {
					continue;
				} 
				if (mark(randomNbr, graph, true)) {
					break;
				}
			}

			cout << "Completed test 1 for N = " << N << endl;
			assertCompletion(1, graph);
			resetGraph(graph);

			//Sends random node not sent before
			random_shuffle(v2.begin(), v2.end());
			while (!v2.empty()) {
				int i = v2.back();
				v2.pop_back();
				++itr2;
				if (!graph[i].marked) {
					if (mark(i, graph, true)) {
						break;
					}
				}
			}

			cout << "Completed test 2 for N = " << N << endl;
			assertCompletion(2, graph);
			resetGraph(graph);

			random_shuffle(v3.begin(), v3.end());
			for (int i : v3) {
				if (!graph[i].marked) {
					++itr3;
					if (mark(i, graph, true)) {
						break;
					}
				}
			}

			cout << "Completed test 3 for N = " << N << endl;
			assertCompletion(3, graph);

			iterations1[r][l - 2] += itr1;
			iterations2[r][l - 2] += itr2;
			iterations3[r][l - 2] += itr3;
		}
	}
	vector<double> means1;
	vector<double> means2;
	vector<double> means3;
	for (int i = 0; i < MAX_LAYERS - 1; ++i) {
		double mean1 = 0.0;
		double mean2 = 0.0;
		double mean3 = 0.0;
		for (int j = 0; j < N_ROUNDS; ++j) {
			mean1 += iterations1[j][i];
			mean2 += iterations2[j][i];
			mean3 += iterations3[j][i];
		}
		mean1 /= N_ROUNDS;
		mean2 /= N_ROUNDS;
		mean3 /= N_ROUNDS;
		means1.push_back(mean1);
		means2.push_back(mean2);
		means3.push_back(mean3);
		double dev1 = 0.0;
		double dev2 = 0.0;
		double dev3 = 0.0;
		for (int j = 0; j < N_ROUNDS; ++j) {
			dev1 += pow(mean1 - iterations1[j][i], 2);
			dev2 += pow(mean2 - iterations2[j][i], 2);
			dev3 += pow(mean3 - iterations3[j][i], 2);
		}
		dev1 /= (pow(2, i + 2) - 1);
		dev1 = sqrt(dev1);
		dev2 /= (pow(2, i + 2) - 1);
		dev2 = sqrt(dev2);
		dev3 /= (pow(2, i + 2) - 1);
		dev3 = sqrt(dev3);
		ofs << "N = " << pow(2, i + 2) - 1 << endl;
		ofs << "Method 1 mean: " << mean1 << ", deviation: " << dev1 << endl;
		ofs << "Method 2 mean: " << mean2 << ", deviation: " << dev2 << endl;
		ofs << "Method 3 mean: " << mean3 << ", deviation: " << dev3 << endl;
		ofs << endl;
	}

	ofs.close();

}