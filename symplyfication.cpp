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




using namespace std;


typedef long double uint64;


const uint32_t FIRST_VERTEX_IN_ADJACENT_LIST = 0;
const bool LOCK_SUITORS_QUEUE = true;
const uint32_t NO_MORE_OFFER = 2000000;
const uint32_t NO_HEAVY_LIGHT_EDGE = 0;


unsigned int bvalue(unsigned int method, unsigned long node_id) {
switch (method) {
default: return (2 * node_id + method) % 10;
case 0: return 4;
case 1: return 7;
} 
}



struct Edge
{
	uint32_t neighbour;
	uint32_t originalNeighbour;
	uint32_t size; 	
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


uint32_t numberOfEdges;
uint32_t liczbaWatkow;
string plikGrafu;	
uint32_t limit_b;

vector<uint32_t> originalNodesIds;
uint32_t numberOfVertices;

queue<uint32_t> Q;
unordered_set<uint32_t> verticesToRepeat;
std::mutex lock_Q;
std::mutex lock_verticesToRepeat;

vector< vector<Edge> > adjList;
vector<uint32_t> posInAdjList;


vector<SuitorsQueue> suitors;
vector<uint32_t> sendProposals;//ZMIANA
vector<uint32_t> bValue;

vector<std::mutex> suitorsProtection;
vector<std::mutex> sendProposalsProtection;


vector<uint32_t> edges;
vector<uint32_t> partialSortUpperBound;
vector<uint32_t> pValue;

bool neighboursExhausted(uint32_t vertex) {
	uint32_t pos = posInAdjList[vertex];
	return pos == adjList[vertex].size();
}


Edge getLowestOffer(uint32_t vertex) {
	if(bValue[vertex] == 0) return Edge{0,0,NO_MORE_OFFER};
	if(suitors[vertex].size() < bValue[vertex]) {
		return Edge{0,0,0};
	}
	return suitors[vertex].top();
}


bool isEligiblePartner(uint32_t u, const Edge& adjPartner, bool lockSuitorsQueue) {
	Edge lowestOffer;
	Edge newOffer = Edge{u, originalNodesIds[u], adjPartner.size};
	if(lockSuitorsQueue) suitorsProtection[adjPartner.neighbour].lock();
	lowestOffer = getLowestOffer(adjPartner.neighbour);
	bool isEligible = edgeGreater(newOffer, lowestOffer);
	if(lockSuitorsQueue) suitorsProtection[adjPartner.neighbour].unlock();
	return isEligible;
}


void sortNextBatch(uint32_t vertex) {
	uint32_t start = partialSortUpperBound[vertex];
	uint32_t end = start + pValue[vertex];
	auto potentialLast = adjList[vertex].begin()+end;
	auto last = (potentialLast < adjList[vertex].end()) ? potentialLast : adjList[vertex].end()-1;
	partial_sort (adjList[vertex].begin()+start, last, adjList[vertex].end(), edgeGreater);
	partialSortUpperBound[vertex] += pValue[vertex];
}


Edge* findEligiblePartner(uint32_t vertex) {
	uint32_t pos = posInAdjList[vertex];
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


uint32_t currentProposalsNumber(uint32_t vertex) {
	uint32_t proposals;
	sendProposalsProtection[vertex].lock();
	proposals = sendProposals[vertex];
	sendProposalsProtection[vertex].unlock();
	return proposals;
}


bool addToSuitorsQueue(uint32_t& annuledVertex,const uint32_t& u,const uint32_t& edgeSize,const uint32_t& p) {
	Edge e{u, originalNodesIds[u], edgeSize};
	bool proposalAnnuled = false;

	if(suitors[p].size() >= bValue[p]) {
		annuledVertex = suitors[p].top().neighbour;
		suitors[p].pop();
		proposalAnnuled = true;
	}
	suitors[p].push(e);
	return proposalAnnuled;
}


bool makeUsuitorOfP(uint32_t& annuledVertex, uint32_t u, uint32_t edgeSize, uint32_t p) {
	return addToSuitorsQueue(annuledVertex, u, edgeSize, p);
}


void updateAnnuledVertex(const bool& proposalAnnuled, const int& annuledVertex, const int& vertexNoLongerSuitet) {
	if(!proposalAnnuled) return;
	
	lock_verticesToRepeat.lock();
	verticesToRepeat.insert(annuledVertex);
	lock_verticesToRepeat.unlock();

	sendProposalsProtection[annuledVertex].lock();
	sendProposals[annuledVertex]--;
	sendProposalsProtection[annuledVertex].unlock();
}


void addToProposals(uint32_t vertex, uint32_t partner) {
	sendProposalsProtection[vertex].lock();
	sendProposals[vertex]++;
	sendProposalsProtection[vertex].unlock();
}


void makeProposes(uint32_t vertex) {
	uint32_t proposes = currentProposalsNumber(vertex); //uzywaj size mapy
	Edge* partner;
	uint32_t annuledVertex;
	bool proposalAnnuled;
	bool createdProposal;
	
	while (proposes < bValue[vertex] && !neighboursExhausted(vertex)) {
		createdProposal = false;
		partner = findEligiblePartner(vertex);
		if (partner != EXAUSTED) {
			suitorsProtection[partner->neighbour].lock();
			if(isEligiblePartner(vertex, *partner, !LOCK_SUITORS_QUEUE)) {
				proposalAnnuled = makeUsuitorOfP(annuledVertex, vertex, partner->size, partner->neighbour);
				createdProposal = true;
				proposes++;
			}
			suitorsProtection[partner->neighbour].unlock();
			if (createdProposal) {
				addToProposals(vertex, partner->neighbour);
				updateAnnuledVertex(proposalAnnuled, annuledVertex, partner->neighbour);
			}
		}		
	}
}

bool nextPartnerIsHeavy(uint32_t v, uint32_t heavyEdgePivot) {
	uint32_t pos = posInAdjList[v];
	return adjList[v][pos].size >= heavyEdgePivot;
}


void makeProposesFirstTime(uint32_t vertex, uint32_t heavyEdgePivot) {
	uint32_t proposes = currentProposalsNumber(vertex); 
	Edge* partner;
	uint32_t annuledVertex;
	bool createdProposal;
	uint32_t pos;
	bool proposalAnnuled;
	
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
				proposalAnnuled = makeUsuitorOfP(annuledVertex, vertex, partner->size, partner->neighbour);
				createdProposal = true;
				proposes++;
			}
			suitorsProtection[partner->neighbour].unlock();
			if (createdProposal) {
				addToProposals(vertex, partner->neighbour);
				updateAnnuledVertex(proposalAnnuled, annuledVertex, partner->neighbour);
			}
		}		
	}
}



uint32_t getNextVertex(bool& haveNextVertex) {
	uint32_t vertex = 0;
	haveNextVertex = false;
	lock_Q.lock();
	if(!Q.empty()) {
		vertex = Q.front();
		Q.pop();
		haveNextVertex = true;
	}
	lock_Q.unlock();
	return vertex;	
}

void processVertex(uint32_t vertex) {
	bool haveNextVertex = true;
	while(haveNextVertex) {
		makeProposes(vertex);
		vertex = getNextVertex(haveNextVertex);
	}
}

void processVertexFirstTime(uint32_t vertex, uint32_t heavyEdgePivot) {
	bool haveNextVertex = true;
	while(haveNextVertex) {
		makeProposesFirstTime(vertex, heavyEdgePivot);
		vertex = getNextVertex(haveNextVertex);
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



void process_Q(uint32_t heavyEdgePivot) {
	vector<std::thread> threads;
	uint32_t nextVertex;
	bool haveNextVertex;
	for(uint32_t i = 1; i <= liczbaWatkow; i++) {
		nextVertex = getNextVertex(haveNextVertex);
		if (!haveNextVertex) break;
		if(i == liczbaWatkow) {
			if(heavyEdgePivot > 0) processVertexFirstTime(nextVertex, heavyEdgePivot);
			else processVertex(nextVertex);
		}
		else {
			if(heavyEdgePivot > 0) threads.emplace_back(thread(processVertexFirstTime, nextVertex, heavyEdgePivot));
			else threads.emplace_back(thread(processVertex, nextVertex));
		}
	}
	join_all(threads);
	for (const auto& vertex: verticesToRepeat)
   		Q.push(vertex);
	verticesToRepeat.clear();
}




void setBPvalues(uint32_t generator) {
	static uint32_t batch_size_parameter = 9;
	for(uint32_t i = 0; i < bValue.size(); i++) {
		bValue[i] = bvalue(generator, originalNodesIds[i]);
		pValue[i] = batch_size_parameter*bValue[i];
	}	
}



void clearSuitorsAssociatedSetsAndQueues() {
	for (SuitorsQueue& queue: suitors)
   		while(!queue.empty()) queue.pop();

	for(uint32_t i = 0; i < posInAdjList.size(); i++) {
		sendProposals[i] = 0;
		posInAdjList[i] = FIRST_VERTEX_IN_ADJACENT_LIST;
	}
}


void moveAllVerticesTo_Q() {
	while(!Q.empty()) Q.pop();
	//heura
	for(int i = adjList.size()-1; i >= 0; i--)
		Q.push((uint32_t)i);	
}


uint32_t sumOfvertexSuitors(uint32_t vertex, uint32_t &ilosc) {
	uint32_t sum = 0;

	while(!suitors[vertex].empty()) {
		sum += suitors[vertex].top().size;
		suitors[vertex].pop();
		ilosc++;
	}
	return sum;
}


uint32_t sumOfAllSuitorsEdges() {
	uint32_t sum = 0;
	uint32_t ilosc = 0;

	for(uint32_t vertex = 0; vertex < suitors.size(); vertex++) {
		sum += sumOfvertexSuitors(vertex, ilosc);
	}
	return sum;
}

uint32_t findValueOfbMatching(uint32_t numberOfVertices, uint32_t generator, uint32_t heavyEdgesPivotVal) {
	setBPvalues(generator);
	clearSuitorsAssociatedSetsAndQueues();
	moveAllVerticesTo_Q();	
	process_Q(heavyEdgesPivotVal);
	while(!Q.empty()) {
		process_Q(0);
	}
	return sumOfAllSuitorsEdges();
}


uint32_t findEdgeAverage() {
	float average = 0;
	uint64 sum = 0;
	uint32_t counter = 0;
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
	return (uint32_t)(average+1);
}



uint32_t findUniversalEdgePivot() {
	uint32_t k = 5;
	uint32_t d = 8;
	uint32_t c = 2;
	uint32_t B = adjList.size()*k;
	auto universalPivot = edges.begin()+B*k;
	nth_element(edges.begin(), universalPivot, edges.end());
	while(k > 0) {
		uint32_t dif = B;
		if(universalPivot < edges.end()) {
			nth_element (edges.begin(), universalPivot, edges.end());
			return *(universalPivot)+findEdgeAverage()/d+c;	
		}
		k--;
		universalPivot = edges.begin()+B*k;
	}	
	return NO_HEAVY_LIGHT_EDGE;
}


void addEdgeCount(uint32_t readVertex, map<uint32_t,uint32_t>& renamedVertex, vector<uint32_t>& vertexAdjacentEdges) {
	auto it = renamedVertex.find(readVertex);
	if (it != renamedVertex.end()) {
		vertexAdjacentEdges[it->second]++;
	}
	else {
		renamedVertex.insert ( pair<uint32_t,uint32_t>(readVertex,renamedVertex.size()) );
		originalNodesIds.emplace_back(readVertex);
		vertexAdjacentEdges.emplace_back(1);
	}
}


void getValuesFromLine(uint32_t& v1, uint32_t& v2, uint32_t& size, string& line) {
	string tmp = "";
	uint32_t which = 0;
	for(uint32_t i = 0; i < line.length(); i++) {
		if(line[i] == ' ') {
			if (which == 0) {v1 = atoi(tmp.c_str()); which++;}
			else v2 = atoi(tmp.c_str());
			tmp = "";
		}
		else tmp += line[i];	
	}
	size = atoi(tmp.c_str());
}


void renameAndCountEdges(string fileName, map<uint32_t,uint32_t>& renamedVertex, vector<uint32_t>& vertexAdjacentEdges) {
	string line;
	numberOfEdges = 0;

	ifstream graphFile (fileName);
	uint32_t v1, v2, size;
	if (graphFile.is_open())
	{
		while ( getline (graphFile,line) )
		{
			if(line[0] == '#') continue;
			getValuesFromLine(v1, v2, size, line);
			numberOfEdges++;
			addEdgeCount(v1, renamedVertex, vertexAdjacentEdges);
			addEdgeCount(v2, renamedVertex, vertexAdjacentEdges);			
		}
		graphFile.close();
	}
}


void allocateAdjacentLists(vector<uint32_t>& vertexAdjacentEdges) {
	adjList.resize(vertexAdjacentEdges.size(), vector<Edge> (0));
	for(uint32_t i = 0; i < vertexAdjacentEdges.size(); i++) {
		adjList[i].resize(vertexAdjacentEdges[i], Edge {0,0,0});
	}
}


void addEdgeToAdjacentList(uint32_t vertex, uint32_t neighbour, uint32_t size, vector<uint32_t>& vertexAdjacentEdges) {
	uint32_t pos = vertexAdjacentEdges[vertex];
	vertexAdjacentEdges[vertex]--;
	adjList[vertex][pos-1].neighbour = neighbour;
	adjList[vertex][pos-1].originalNeighbour = originalNodesIds[neighbour];	
	adjList[vertex][pos-1].size = size;
}


void createAdjacentLists(string fileName, map<uint32_t,uint32_t>& renamedVertex, vector<uint32_t>& vertexAdjacentEdges, vector<uint32_t>& edges) {
	string line;
	ifstream graphFile (fileName);
	uint32_t v1, v2, size;
	uint32_t renamed_v1, renamed_v2;
	uint32_t edgeCount = 0;	
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
	uint32_t numberOfVertices = adjList.size();
	bValue.resize(numberOfVertices, 0);
	suitors.resize(numberOfVertices);
	posInAdjList.resize(numberOfVertices, 0);
	sendProposals.resize(numberOfVertices, 0);
	partialSortUpperBound.resize(numberOfVertices, 0);
	pValue.resize(numberOfVertices, 0);
}



int main(int argc, char **argv)
{
	map<uint32_t,uint32_t> renamedVertex;
	vector<uint32_t> vertexAdjacentEdges;
	
	liczbaWatkow = atoi(argv[1]);
	plikGrafu = argv[2];	
	limit_b = atoi(argv[3]);

	renameAndCountEdges(plikGrafu, renamedVertex, vertexAdjacentEdges);
	createAdjacentLists(plikGrafu, renamedVertex, vertexAdjacentEdges, edges);
	allocateProtection();
	allocateControlConteners();
	
	uint32_t maximazedSum;
	uint32_t universalEdgePivot = findUniversalEdgePivot();	
	
	for(uint32_t generator = 0; generator <= limit_b; generator++) {
		maximazedSum = findValueOfbMatching(numberOfVertices, generator, universalEdgePivot)/2;
		cout<<maximazedSum<<endl;
	}

    return 0;                                                                  
}
