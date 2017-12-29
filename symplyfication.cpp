#include <iostream>
#include <thread>
#include <chrono>
#include <vector>
#include <atomic>
#include <mutex>
#include <queue>
#include <unordered_set>
#include <algorithm>
#include <map>
#include <utility>
#include <fstream>
#include <numeric>
#include <limits>

#define EXAUSTED nullptr
#define NOT_ANNULED -1



using namespace std;



typedef long double uint64;
typedef unsigned int uint;

const int LACK_OF_VERTICES = -1;
const int FIRST_VERTEX_IN_ADJACENT_LIST = 0;
const bool LOCK_SUITORS_QUEUE = true;
const int NO_MORE_OFFER = 2000000;
const int NO_HEAVY_LIGHT_EDGE = 0;




unsigned int bvalue(unsigned int method, unsigned long node_id) {
switch (method) {
default: return (2 * node_id + method) % 10;
case 0: return 4;
case 1: return 7;
} 
}




struct Edge
{
	int neighbour;
	int originalNeighbour;
	int size; 	
};


class edgeComparision
{
  bool reverse;
public:
  edgeComparision(const bool& revparam=true)
    {reverse=revparam;}

  bool operator() (const Edge& lhs, const Edge& rhs) const
  {
	if (lhs.size == rhs.size) {
		if(reverse) return lhs.originalNeighbour > rhs.originalNeighbour;
		return lhs.originalNeighbour < rhs.originalNeighbour;	
	}
	if(reverse) {
		return lhs.size > rhs.size;
	}
		return lhs.size < rhs.size;
  	}
};

typedef priority_queue<Edge,vector<Edge>,edgeComparision> SuitorsQueue;


bool edgeGreater(const Edge& e1, const Edge& e2) {//TO DO porownanie starych
	if (e1.size > e2.size) return true;
	else if(e1.size == e2.size) return e1.originalNeighbour > e2.originalNeighbour;
	return false;
}


int numberOfEdges;
int liczbaWatkow;
string plikGrafu;	
int limit_b;

vector<int> originalNodesIds;
int numberOfVertices;

queue<int> Q;
unordered_set<int> verticesToRepeat;
std::mutex lock_Q;
std::mutex lock_verticesToRepeat;

vector< vector<Edge> > adjList;
vector<int> posInAdjList;


vector<SuitorsQueue> suitors;
//vector< unordered_set<int> > sendProposals;
vector<int> sendProposals;
vector<int> bValue;

vector<std::mutex> suitorsProtection;
vector<std::mutex> sendProposalsProtection;


vector<int> edges;
vector<int> partialSortUpperBound;
vector<int> pValue;

bool neighboursExhausted(int vertex) {
	int pos = posInAdjList[vertex];
	return pos == adjList[vertex].size();
}


Edge getLowestOffer(int vertex) {
	if(bValue[vertex] == 0) return Edge{0,0,NO_MORE_OFFER};
	if(suitors[vertex].size() < bValue[vertex]) {
		return Edge{0,0,0};
	}
	return suitors[vertex].top();
}


bool isEligiblePartner(int u, const Edge& adjPartner, bool lockSuitorsQueue) {
	Edge lowestOffer;
	Edge newOffer = Edge{u, originalNodesIds[u], adjPartner.size};
	if(lockSuitorsQueue) suitorsProtection[adjPartner.neighbour].lock();
	lowestOffer = getLowestOffer(adjPartner.neighbour);
	bool isEligible = edgeGreater(newOffer, lowestOffer);
	if(lockSuitorsQueue) suitorsProtection[adjPartner.neighbour].unlock();
	return isEligible;
}


void sortNextBatch(int vertex) {
	int start = partialSortUpperBound[vertex];
	int end = start + pValue[vertex];
	auto potentialLast = adjList[vertex].begin()+end;
	auto last = (potentialLast < adjList[vertex].end()) ? potentialLast : adjList[vertex].end()-1;
	partial_sort (adjList[vertex].begin()+start, last, adjList[vertex].end(), edgeGreater);
	partialSortUpperBound[vertex] += pValue[vertex];
}


Edge* findEligiblePartner(int vertex) {
	int pos = posInAdjList[vertex];
	while (pos < adjList[vertex].size()) {
		if(pos >= partialSortUpperBound[vertex]) sortNextBatch(vertex);
		if(isEligiblePartner(vertex, adjList[vertex][pos], LOCK_SUITORS_QUEUE)) {
			posInAdjList[vertex]++;			
			return &adjList[vertex][pos];
		}
		pos++;
		posInAdjList[vertex]++;	
	}
	return EXAUSTED;
}


int currentProposalsNumber(int vertex) {
	int proposals;
	sendProposalsProtection[vertex].lock();
	proposals = sendProposals[vertex];
	sendProposalsProtection[vertex].unlock();
	return proposals;
}


int addToSuitorsQueue(int u, int edgeSize, int p) {
	Edge e{u, originalNodesIds[u], edgeSize};
	int annuledVertex = NOT_ANNULED;

	if(suitors[p].size() >= bValue[p]) {
		annuledVertex = suitors[p].top().neighbour;
		suitors[p].pop();
	}
	suitors[p].push(e);
	return annuledVertex;
}


int makeUsuitorOfP(int u, int edgeSize, int p) {
	int annuledVertex = addToSuitorsQueue(u, edgeSize, p);
	return annuledVertex;
}


void updateAnnuledVertex(int annuledVertex, int vertexNoLongerSuitet) {
	if(annuledVertex == NOT_ANNULED) return;
	
	lock_verticesToRepeat.lock();
	verticesToRepeat.insert(annuledVertex);
	lock_verticesToRepeat.unlock();

	sendProposalsProtection[annuledVertex].lock();
	sendProposals[annuledVertex]--;
	sendProposalsProtection[annuledVertex].unlock();
}


void addToProposals(int vertex, int partner) {
	sendProposalsProtection[vertex].lock();
	sendProposals[vertex]++;
	sendProposalsProtection[vertex].unlock();
}


void makeProposes(int vertex) {
	int proposes = currentProposalsNumber(vertex);
	Edge* partner;
	int annuledVertex;
	bool createdProposal;
	
	while (proposes < bValue[vertex] && !neighboursExhausted(vertex)) {
		createdProposal = false;
		partner = findEligiblePartner(vertex);
		if (partner != EXAUSTED) {
			suitorsProtection[partner->neighbour].lock();
			if(isEligiblePartner(vertex, *partner, !LOCK_SUITORS_QUEUE)) {
				annuledVertex = makeUsuitorOfP(vertex, partner->size, partner->neighbour);
				createdProposal = true;
				proposes++;
			}
			if (createdProposal) {
				//addToProposals(vertex, partner->neighbour);
			}
			suitorsProtection[partner->neighbour].unlock();
			if (createdProposal) {
				addToProposals(vertex, partner->neighbour);
				updateAnnuledVertex(annuledVertex, partner->neighbour);
			}
		}		
	}
}

bool nextPartnerIsHeavy(int v, int heavyEdgePivot) {
	int pos = posInAdjList[v];
	return adjList[v][pos].size >= heavyEdgePivot;
}


void makeProposesFirstTime(int vertex, int heavyEdgePivot) {
	int proposes = currentProposalsNumber(vertex);
	Edge* partner;
	int annuledVertex;
	bool createdProposal;
	int pos;
	
	while (proposes < bValue[vertex] && !neighboursExhausted(vertex)) {
		
		if(!nextPartnerIsHeavy(vertex, heavyEdgePivot)) {
			lock_verticesToRepeat.lock();
			verticesToRepeat.insert(vertex);
			lock_verticesToRepeat.unlock();
			break;
		}

		createdProposal = false;
		partner = findEligiblePartner(vertex);
		if (partner != EXAUSTED) {
			suitorsProtection[partner->neighbour].lock();
			if(isEligiblePartner(vertex, *partner, !LOCK_SUITORS_QUEUE)) {
				annuledVertex = makeUsuitorOfP(vertex, partner->size, partner->neighbour);
				createdProposal = true;
				proposes++;
			}
			if (createdProposal) {
				//addToProposals(vertex, partner->neighbour);
			}
			suitorsProtection[partner->neighbour].unlock();
			if (createdProposal) {
				addToProposals(vertex, partner->neighbour);
				updateAnnuledVertex(annuledVertex, partner->neighbour);
			}
		}		
	}
}



int getNextVertex() {
	int vertex;
	lock_Q.lock();
		if(!Q.empty()) {
			vertex = Q.front();
			Q.pop();
		}
		else vertex = LACK_OF_VERTICES;
		lock_Q.unlock();
	return vertex;	
}

void processVertex(int vertex) {
	while(vertex != LACK_OF_VERTICES) {
		makeProposes(vertex);
		vertex = getNextVertex();
	}
}

void processVertexFirstTime(int vertex, int heavyEdgePivot) {
	while(vertex != LACK_OF_VERTICES) {
		makeProposesFirstTime(vertex, heavyEdgePivot);
		vertex = getNextVertex();
	}
}


void do_join(std::thread& t)
{
    t.join();
}


void join_all(vector<std::thread>& v)
{
    for_each(v.begin(),v.end(),do_join);
}





void process_Q(int heavyEdgePivot) {
	vector<std::thread> threads;
	int nextVertex;
	for(int i = 1; i <= liczbaWatkow; i++) {
		nextVertex = getNextVertex();
		if (nextVertex == LACK_OF_VERTICES) break;
		if(i == liczbaWatkow) {
			if(heavyEdgePivot > 0) processVertexFirstTime(nextVertex, heavyEdgePivot);
			else processVertex(nextVertex);//process vertex yourself
		}
		else {
			if(heavyEdgePivot > 0) threads.emplace_back(thread(processVertexFirstTime, nextVertex, heavyEdgePivot));
			else threads.emplace_back(thread(processVertex, nextVertex));//create process to do the job
		}
	}
	join_all(threads);
	for (const auto& vertex: verticesToRepeat)
   		Q.push(vertex);
	verticesToRepeat.clear();
}




void setBPvalues(int generator, vector<int>& originalNodesIds) {
	static const int batch_size_parameter = 9;
	for(int i = 0; i < bValue.size(); i++) {
		bValue[i] = bvalue(generator, originalNodesIds[i]);
		pValue[i] = batch_size_parameter*bValue[i];
	}	
	for(int i = 0; i < bValue.size(); i++) {
		//cout<<"oryginalnie "<<originalNodesIds[i]<<" "<<bValue[i]<<endl;
	}	
}



void clearSuitorsAssociatedSetsAndQueues() {
	for (SuitorsQueue& queue: suitors)
   		while(!queue.empty()) queue.pop();
	for(vector<int>::iterator it = sendProposals.begin(); it != sendProposals.end(); ++it) {
		*(it) = 0;
	}	
	for(int i = 0; i < posInAdjList.size(); i++) 
		posInAdjList[i] = FIRST_VERTEX_IN_ADJACENT_LIST;
	}


void moveAllVerticesTo_Q() {
	while(!Q.empty()) Q.pop();
	//for(int i = 0; i < adjList.size(); i++) 
	//	Q.push(i);	

	//heura
	for(int i = adjList.size()-1; i >= 0; i--)
		Q.push(i);	
}


int sumOfvertexSuitors(int vertex, int &ilosc) {
	int sum = 0;

	while(!suitors[vertex].empty()) {
		sum += suitors[vertex].top().size;
		suitors[vertex].pop();
		ilosc++;
	}
	return sum;
}


int sumOfAllSuitorsEdges() {
	int sum = 0;
	int ilosc = 0;

	for(int vertex = 0; vertex < suitors.size(); vertex++) {
		sum += sumOfvertexSuitors(vertex, ilosc);
	}
	return sum;
}

int findValueOfbMatching(int numberOfVertices, int generator, vector<int>& originalNodesIds, int heavyEdgesPivotVal) {
	setBPvalues(generator, originalNodesIds);
	//cout<<"heavyEdgePivot: "<<heavyEdgesPivotVal<<endl;
	clearSuitorsAssociatedSetsAndQueues();
	if(!Q.empty())
		cout<<"JA PIERDOLE\n";
	
	moveAllVerticesTo_Q();
	
	process_Q(heavyEdgesPivotVal);
	while(!Q.empty()) {
		process_Q(0);
	}
	return sumOfAllSuitorsEdges();
}


int findEdgeAverage() {
	float average = 0;
	uint64 sum = 0;
	int counter = 0;
	uint64 limit = 10000000;
	
	for (auto& size: edges) {
		if(counter > limit) {
			average += (float) (sum)/counter;
			counter = 0;
		}
   		sum += size;
		counter++;
	}
	average += (float) (sum)/counter;
	return (int)(average+1);
}


int heavyEdgesPivotVal(int k) {
	uint64 init = 0;
	uint64 B = accumulate(bValue.begin(),bValue.end(),init);
	while(k > 0) {
		uint64 dif = B;
		if(edges.begin()+k*B < edges.end()) {
			nth_element (edges.begin(), edges.begin()+k*B, edges.end());
			return *(edges.begin()+k*B);	
		}
		k--;
	}
	return NO_HEAVY_LIGHT_EDGE;
}

int findUniversalEdgePivot() {
	//oblicz dajac kazdemu wierzcholkowi 3 krawedzie
	int k = 5;
	int d = 8;
	int c = 2;
	int B = adjList.size()*k;
	cout<<"B "<<B<<endl;
	auto universalPivot = edges.begin()+B*k;
	nth_element(edges.begin(), universalPivot, edges.end());
	while(k > 0) {
		int dif = B;
		if(universalPivot < edges.end()) {
			nth_element (edges.begin(), universalPivot, edges.end());
			return *(universalPivot)+findEdgeAverage()/d+c;	
		}
		k--;
		universalPivot = edges.begin()+B*k;
	}	
	return NO_HEAVY_LIGHT_EDGE;
}


void addEdgeCount(int readVertex, map<int,int>& renamedVertex, vector<int>& vertexAdjacentEdges) {
	auto it = renamedVertex.find(readVertex);
	if (it != renamedVertex.end()) {
		vertexAdjacentEdges[it->second]++;
	}
	else {
		renamedVertex.insert ( pair<int,int>(readVertex,renamedVertex.size()) );
		originalNodesIds.emplace_back(readVertex);
		vertexAdjacentEdges.emplace_back(1);
	}
}


void getValuesFromLine(int& v1, int& v2, int& size, string& line) {
	string tmp = "";
	int which = 0;
	for(int i = 0; i < line.length(); i++) {
		if(line[i] == ' ') {
			if (which == 0) {v1 = atoi(tmp.c_str()); which++;}
			else v2 = atoi(tmp.c_str());
			tmp = "";
		}
		else tmp += line[i];	
	}
	size = atoi(tmp.c_str());
}


void renameAndCountEdges(string fileName, map<int,int>& renamedVertex, vector<int>& vertexAdjacentEdges) {
	string line;
	numberOfEdges = 0;

	ifstream graphFile (fileName);
	int v1, v2, size;
	if (graphFile.is_open())
	{
		while ( getline (graphFile,line) )
		{

			//cout<<"uzyskana linia "<<line<<"\n";
			if(line[0] == '#') continue;
			getValuesFromLine(v1, v2, size, line);
			numberOfEdges++;
			addEdgeCount(v1, renamedVertex, vertexAdjacentEdges);
			addEdgeCount(v2, renamedVertex, vertexAdjacentEdges);			
		}
		graphFile.close();
	}
}


void allocateAdjacentLists(vector<int>& vertexAdjacentEdges) {
	adjList.resize(vertexAdjacentEdges.size(), vector<Edge> (0));
	for(int i = 0; i < vertexAdjacentEdges.size(); i++) {
		adjList[i].resize(vertexAdjacentEdges[i]);
	}
}


void addEdgeToAdjacentList(int vertex, int neighbour, int size, vector<int>& vertexAdjacentEdges) {
	int pos = vertexAdjacentEdges[vertex];
	vertexAdjacentEdges[vertex]--;
	adjList[vertex][pos-1].neighbour = neighbour;
	adjList[vertex][pos-1].originalNeighbour = originalNodesIds[neighbour];	
	adjList[vertex][pos-1].size = size;
}


void createAdjacentLists(string fileName, map<int,int>& renamedVertex, vector<int>& vertexAdjacentEdges, vector<int>& edges) {
	string line;

	ifstream graphFile (fileName);
	int v1, v2, size;
	int renamed_v1, renamed_v2;
	int edgeCount = 0;	
	cout<<"number of edges"<<numberOfEdges<<endl;
	edges.resize(numberOfEdges);
	allocateAdjacentLists(vertexAdjacentEdges);
	
	if (graphFile.is_open())
	{
		while ( getline (graphFile,line) )
		{
			if(line[0] == '#') continue;
			getValuesFromLine(v1, v2, size, line);
			renamed_v1 = renamedVertex[v1];
			renamed_v2 = renamedVertex[v2];
			edges[edgeCount++] = size;
			
			addEdgeToAdjacentList(renamed_v1, renamed_v2, size, vertexAdjacentEdges);
			addEdgeToAdjacentList(renamed_v2, renamed_v1, size, vertexAdjacentEdges);			
		}
		graphFile.close();
	}
	vertexAdjacentEdges.clear();
}


void allocateProtection() {	
	suitorsProtection = vector<std::mutex> (adjList.size());
	sendProposalsProtection = vector<std::mutex> (adjList.size());
}


void allocateControlConteners() {
	int numberOfVertices = adjList.size();
	bValue.resize(numberOfVertices, 0);
	suitors.resize(numberOfVertices);
	posInAdjList.resize(numberOfVertices, 0);
	sendProposals.resize(numberOfVertices, 0);
	partialSortUpperBound.resize(numberOfVertices, 0);
	pValue.resize(numberOfVertices, 0);
}


void tymczasowy_sort_sekwencyjny() {
	for(int i= 0; i < adjList.size(); i++) {
			std::sort (adjList[i].begin(), adjList[i].end(), edgeGreater);
	}	
}


void drukujAdjacentList() {
	//cout<<"\n\n";
	//cout<<"{neighbour, size}\n";
	for(int i = 0; i < adjList.size(); i++) {
		cout<<"Lista "<<i<<": ";
		for(int j = 0; j < adjList[i].size(); j++) {
			cout<<"{"<<adjList[i][j].originalNeighbour<<" , "<<adjList[i][j].size<<"}"<<"  ";
		}
		cout<<endl;
	}
	//cout<<"wydrukowalem\n";
}






int main(int argc, char **argv)
{
	map<int,int> renamedVertex;
	vector<int> vertexAdjacentEdges;
	
	liczbaWatkow = atoi(argv[1]);
	plikGrafu = argv[2];	
	limit_b = atoi(argv[3]);

	cout<<"zaczynam wczytywac i alokowac \n";
	renameAndCountEdges(plikGrafu, renamedVertex, vertexAdjacentEdges);//tu? rename
	createAdjacentLists(plikGrafu, renamedVertex, vertexAdjacentEdges, edges);//rename poprawny tu	
	allocateProtection();
	cout<<"koncze czytac, zaczynam alokowac \n";
	allocateControlConteners();
	cout<<"konczy wczytywac i alokowac \n";

	
	int maximazedSum;
	int universalEdgePivot = findUniversalEdgePivot();
	cout<<"universalEdgePivot: "<<universalEdgePivot<<endl;
	for(int generator = 0; generator <= limit_b; generator++) {
		maximazedSum = findValueOfbMatching(numberOfVertices, generator, originalNodesIds, universalEdgePivot)/2;
		cout<<maximazedSum<<endl;
	}
	

    return 0;                                                                  
}
