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

#define EXAUSTED nullptr
#define NOT_ANNULED -1



using namespace std;

/*
ULEPSZENIA
-prawdopodobnie nie trzeba trzymac w kazdym wierzcholku zbioru do kogo wyslal propozycje, wystarczy tylko ich liczba
-HEAVY-LIGHT edge

*/




const int LACK_OF_VERTICES = -1;
const int FIRST_VERTEX_IN_ADJACENT_LIST = 0;
const bool LOCK_SUITORS_QUEUE = true;
const int NO_MORE_OFFER = 2000000;

/*
unsigned int bvalue(unsigned int method, unsigned long node_id) {

    switch (method) {

    case 0: return 1;

    default: switch (node_id) {

        case 0: return 2;

        case 1: return 2;

        default: return 1;

        }

    }

}

unsigned int bvalue(unsigned int method, unsigned long node_id) {
	return node_id % 2;
}

*/



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
vector< unordered_set<int> > sendProposals;
vector<int> bValue;

vector<std::mutex> suitorsProtection;
vector<std::mutex> sendProposalsProtection;


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
	//cout<<"lowestOffer originalNei: "<<lowestOffer.originalNeighbour<<endl;
	bool isEligible = edgeGreater(newOffer, lowestOffer);
	if(lockSuitorsQueue) suitorsProtection[adjPartner.neighbour].unlock();
	return isEligible;
}


Edge* findEligiblePartner(int vertex) {
	int pos = posInAdjList[vertex];
	while (pos < adjList[vertex].size()) {
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
	proposals = sendProposals[vertex].size();
	sendProposalsProtection[vertex].unlock();
	return proposals;
}


int addToSuitorsQueue(int u, int edgeSize, int p) {
	//cout<<"wkladam original "<<originalNodesIds[u]<<endl;
	Edge e{u, originalNodesIds[u], edgeSize};
	int annuledVertex = NOT_ANNULED;

	if(suitors[p].size() >= bValue[p]) {
		annuledVertex = suitors[p].top().neighbour;
		suitors[p].pop();
	}
	else {
		//cout<<originalNodesIds[u]<<" become suitor of "<<originalNodesIds[p]<<" edgeSize "<<edgeSize<<endl;
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
	int s1, s2;
	s1 = sendProposals[annuledVertex].size();
	sendProposals[annuledVertex].erase(vertexNoLongerSuitet);
	s2 = sendProposals[annuledVertex].size();
	//if(s1 == s2) cout<<"NOE ZESZ KOWA !!! \n";
	sendProposalsProtection[annuledVertex].unlock();
}


void addToProposals(int vertex, int partner) {
	sendProposalsProtection[vertex].lock();
	sendProposals[vertex].insert(partner);
	sendProposalsProtection[vertex].unlock();
}


void makeProposes(int vertex) {
	int proposes = currentProposalsNumber(vertex); //uzywaj size mapy
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
			suitorsProtection[partner->neighbour].unlock();
			if (createdProposal) {
				updateAnnuledVertex(annuledVertex, partner->neighbour);
				addToProposals(vertex, partner->neighbour);
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
	//bool haveNextVertex = true;
	while(vertex != LACK_OF_VERTICES) {
		makeProposes(vertex);
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



void process_Q() {
	vector<std::thread> threads;
	int nextVertex;
	for(int i = 1; i <= liczbaWatkow; i++) {
		nextVertex = getNextVertex();
		if (nextVertex == LACK_OF_VERTICES) break;
		if(i == liczbaWatkow) {
			processVertex(nextVertex);//process vertex yourself
		}
		else {
			threads.emplace_back(std::thread(processVertex, nextVertex));//create process to do the job
		}
	}
	join_all(threads);
	for (const auto& vertex: verticesToRepeat)
   		Q.push(vertex);
	verticesToRepeat.clear();
}


void setBValue(int generator, vector<int>& originalNodesIds) {
	for(int i = 0; i < bValue.size(); i++) {
		bValue[i] = bvalue(generator, originalNodesIds[i]);
	}	
	for(int i = 0; i < bValue.size(); i++) {
		//cout<<"oryginalnie "<<originalNodesIds[i]<<" "<<bValue[i]<<endl;
	}	
}



void clearSuitorsAssociatedSetsAndQueues() {
	for (SuitorsQueue& queue: suitors)
   		while(!queue.empty()) queue.pop();

	for (auto& u_set: sendProposals)
   		u_set.clear();	
	
	for(int i = 0; i < posInAdjList.size(); i++) 
		posInAdjList[i] = FIRST_VERTEX_IN_ADJACENT_LIST;
}


void moveAllVerticesTo_Q() {
	while(!Q.empty()) Q.pop();
	for(int i = 0; i < adjList.size(); i++) 
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
	//cout<<"ILOSC "<<ilosc<<endl;
	return sum;
}

int findValueOfbMatching(int numberOfVertices, int generator, vector<int>& originalNodesIds) {
	setBValue(generator, originalNodesIds);
	clearSuitorsAssociatedSetsAndQueues();
	if(!Q.empty())
		cout<<"JA PIERDOLE\n";
	
	moveAllVerticesTo_Q();

	while(!Q.empty()) {
		process_Q();
	}
	return sumOfAllSuitorsEdges();
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
	//cout<<"uzyskane wartosci: "<<v1<<" "<<v2<<" "<<size<<endl;
}


void renameAndCountEdges(string fileName, map<int,int>& renamedVertex, vector<int>& vertexAdjacentEdges) {
	string line;

	ifstream graphFile (fileName);
	int v1, v2, size;
	if (graphFile.is_open())
	{
		while ( getline (graphFile,line) )
		{

			//cout<<"uzyskana linia "<<line<<"\n";
			if(line[0] == '#') continue;
			getValuesFromLine(v1, v2, size, line);
			
			addEdgeCount(v1, renamedVertex, vertexAdjacentEdges);
			addEdgeCount(v2, renamedVertex, vertexAdjacentEdges);			
		}
		graphFile.close();
	}
	
	for(int i = 0; i < vertexAdjacentEdges.size(); i++) {
		//cout<<"ilosc krawedzi "<<i<<" (oryginalnie "<<originalNodesIds[i]<<" ):"<<vertexAdjacentEdges[i]<<endl;
		//cout<<i<<" (oryginalnie "<<originalNodesIds[i]<<endl;		
	}	
}


void allocateAdjacentLists(vector<int>& vertexAdjacentEdges) {
	adjList.resize(vertexAdjacentEdges.size(), vector<Edge> (0));
	for(int i = 0; i < vertexAdjacentEdges.size(); i++) {
		adjList[i].resize(vertexAdjacentEdges[i], Edge {0,0,0});
	}
}


void addEdgeToAdjacentList(int vertex, int neighbour, int size, vector<int>& vertexAdjacentEdges) {
	int pos = vertexAdjacentEdges[vertex];
	vertexAdjacentEdges[vertex]--;
	adjList[vertex][pos-1].neighbour = neighbour;
	adjList[vertex][pos-1].originalNeighbour = originalNodesIds[neighbour];	
	adjList[vertex][pos-1].size = size;
}


void createAdjacentLists(string fileName, map<int,int>& renamedVertex, vector<int>& vertexAdjacentEdges) {
	string line;

	ifstream graphFile (fileName);
	int v1, v2, size;
	int renamed_v1, renamed_v2;

	allocateAdjacentLists(vertexAdjacentEdges);
	
	if (graphFile.is_open())
	{
		while ( getline (graphFile,line) )
		{
			if(line[0] == '#') continue;
			getValuesFromLine(v1, v2, size, line);
			renamed_v1 = renamedVertex[v1];
			renamed_v2 = renamedVertex[v2];
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
	bValue.resize(adjList.size(), 0);
	suitors.resize(adjList.size());
	posInAdjList.resize(adjList.size(), 0);
	sendProposals.resize(adjList.size(), unordered_set<int>());
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

//TO DO!
//poprawic dla b[v] == 0
//wyczyscic plik do debugowania
int main(int argc, char **argv)
{
	map<int,int> renamedVertex;
	vector<int> vertexAdjacentEdges;
	
	liczbaWatkow = atoi(argv[1]);
	plikGrafu = argv[2];	
	limit_b = atoi(argv[3]);

	renameAndCountEdges(plikGrafu, renamedVertex, vertexAdjacentEdges);//tu? rename
	createAdjacentLists(plikGrafu, renamedVertex, vertexAdjacentEdges);//rename poprawny tu	
	allocateProtection();
	allocateControlConteners();


	
	tymczasowy_sort_sekwencyjny();//nie zrzuca pamieci? raczej robia to wierzcholki b[v] == 0

	//drukujAdjacentList();
	//cout<<findValueOfbMatching(1, limit_b, originalNodesIds)/2<<endl;

	int maximazedSum;	
	for(int generator = 0; generator <= limit_b; generator++) {
		maximazedSum = findValueOfbMatching(numberOfVertices, generator, originalNodesIds)/2;
		cout<<maximazedSum<<endl;
	}
	

    return 0;                                                                  
}




