#include <vector>
#include <queue>
#include <algorithm>
#include <limits>
#include <memory>
#include <iostream>
#include <unordered_map>
#include <sstream>
#include <sys/wait.h>

// Forward declarations
class MST;
class MSTAlgorithm;

// Graph class
class Graph
{
public:
    struct Edge
    {
        int from, to, weight;
        Edge(int f, int t, int w) : from(f), to(t), weight(w) {}
    };

    Graph(int vertices) : V(vertices)
    {
        adjacencyList.resize(V);
    }

    Graph() {};

    void addEdge(int from, int to, int weight)
    {
        adjacencyList[from].emplace_back(from, to, weight);
        adjacencyList[to].emplace_back(to, from, weight); // For undirected graph
        this->edgecount++;
    }

    int getVerticesCount() const { return V; }
    int getEdgesCount() const { return edgecount; }
    const std::vector<std::vector<Edge>> &getAdjacencyList() const { return adjacencyList; }

private:
    int V; // Number of vertices
    std::vector<std::vector<Edge>> adjacencyList;
    int edgecount = 0;
};

// MST class
class MST
{
public:
    std::vector<Graph::Edge> mstEdges;
    MST(const Graph &graph) : originalGraph(graph) {}

    void addEdge(int from, int to, int weight)
    {
        mstEdges.emplace_back(from, to, weight);
    }

    int getTotalWeight() const
    {
        int total = 0;
        for (const auto &edge : mstEdges)
        {
            total += edge.weight;
        }
        return total;
    }

    int getLongestDistance() const
    {
        return dfs(0).first;
    }

    int getShortestEdge() const
    {
        if (mstEdges.empty())
            return 0;
        return std::min_element(mstEdges.begin(), mstEdges.end(),
                                [](const Graph::Edge &a, const Graph::Edge &b)
                                {
                                    return a.weight < b.weight;
                                })
            ->weight;
    }

    double getAverageDistance() const
    {
        int V = originalGraph.getVerticesCount();
        std::vector<std::vector<int>> dist(V, std::vector<int>(V, std::numeric_limits<int>::max()));

        // Initialize distances
        for (int i = 0; i < V; ++i)
        {
            dist[i][i] = 0;
            for (const auto &edge : originalGraph.getAdjacencyList()[i])
            {
                dist[i][edge.to] = edge.weight;
            }
        }

        // Floyd-Warshall algorithm
        for (int k = 0; k < V; ++k)
        {
            for (int i = 0; i < V; ++i)
            {
                for (int j = 0; j < V; ++j)
                {
                    if (dist[i][k] != std::numeric_limits<int>::max() &&
                        dist[k][j] != std::numeric_limits<int>::max() &&
                        dist[i][k] + dist[k][j] < dist[i][j])
                    {
                        dist[i][j] = dist[i][k] + dist[k][j];
                    }
                }
            }
        }

        // Calculate average distance
        long long sum = 0;
        int count = 0; // The number of ways to get from vx u to v
        for (int i = 0; i < V; ++i)
        {
            for (int j = i; j < V; ++j)
            {
                if (dist[i][j] != std::numeric_limits<int>::max())
                {
                    sum += dist[i][j];
                    count++;
                }
            }
        }
        return static_cast<double>(sum) / count;
    }

private:
    const Graph &originalGraph;

    std::pair<int, int> dfs(int node, int parent = -1) const {
        int maxDepth = 0;
        int maxNode = node;

        for (const auto& edge : mstEdges) {
            int nextNode = (edge.from == node) ? edge.to : (edge.to == node) ? edge.from : -1;
            if (nextNode != -1 && nextNode != parent) {
                auto [depth, farthestNode] = dfs(nextNode, node);
                depth += edge.weight;
                if (depth > maxDepth) {
                    maxDepth = depth;
                    maxNode = farthestNode;
                }
            }
        }

        return {maxDepth, maxNode};
    }
};

// Abstract base class for MST algorithms
class MSTAlgorithm
{
public:
    virtual std::unique_ptr<MST> findMST(const Graph &graph) = 0;
    virtual ~MSTAlgorithm() = default;
};

// Disjoint Set data structure for Kruskal's and Borůvka's algorithms
class DisjointSet
{
private:
    std::vector<int> parent, rank;

public:
    DisjointSet(int n) : parent(n), rank(n, 0)
    {
        for (int i = 0; i < n; i++)
            parent[i] = i;
    }

    int find(int x)
    {
        if (parent[x] != x)
            parent[x] = find(parent[x]);
        return parent[x];
    }

    void unite(int x, int y)
    {
        int xroot = find(x), yroot = find(y);
        if (xroot == yroot)
            return;
        if (rank[xroot] < rank[yroot])
            parent[xroot] = yroot;
        else if (rank[xroot] > rank[yroot])
            parent[yroot] = xroot;
        else
        {
            parent[yroot] = xroot;
            rank[xroot]++;
        }
    }
};

// Borůvka's Algorithm
class BoruvkaAlgorithm : public MSTAlgorithm
{
public:
    std::unique_ptr<MST> findMST(const Graph &graph) override
    {
        auto mst = std::make_unique<MST>(graph);
        int V = graph.getVerticesCount();
        DisjointSet ds(V);

        // Store the cheapest edge for each component
        std::vector<Graph::Edge> cheapest(V, {-1, -1, std::numeric_limits<int>::max()});

        int numComponents = V;

        while (numComponents > 1)
        {
            // Find the cheapest edge for each component
            for (int i = 0; i < V; ++i)
            {
                for (const auto &edge : graph.getAdjacencyList()[i])
                {
                    int set1 = ds.find(edge.from);
                    int set2 = ds.find(edge.to);

                    if (set1 != set2)
                    {
                        if (edge.weight < cheapest[set1].weight)
                        {
                            cheapest[set1] = edge;
                        }
                        if (edge.weight < cheapest[set2].weight)
                        {
                            cheapest[set2] = edge;
                        }
                    }
                }
            }

            // Add the cheapest edges to the MST
            for (int i = 0; i < V; ++i)
            {
                if (cheapest[i].from != -1)
                {
                    int set1 = ds.find(cheapest[i].from);
                    int set2 = ds.find(cheapest[i].to);

                    if (set1 != set2)
                    {
                        mst->addEdge(cheapest[i].from, cheapest[i].to, cheapest[i].weight);
                        ds.unite(set1, set2);
                        numComponents--;
                    }
                }
            }

            // Reset cheapest edges for next iteration
            std::fill(cheapest.begin(), cheapest.end(), Graph::Edge(-1, -1, std::numeric_limits<int>::max()));
        }

        return mst;
    }
};

// Prim's Algorithm
class PrimAlgorithm : public MSTAlgorithm
{
public:
    std::unique_ptr<MST> findMST(const Graph &graph) override
    {
        auto mst = std::make_unique<MST>(graph);
        int V = graph.getVerticesCount();
        std::vector<int> key(V, std::numeric_limits<int>::max());
        std::vector<bool> inMST(V, false);
        std::vector<int> parent(V, -1);

        std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> pq;

        key[0] = 0;
        pq.push({0, 0});

        while (!pq.empty())
        {
            int u = pq.top().second;
            pq.pop();

            if (inMST[u])
                continue;

            inMST[u] = true;

            for (const auto &edge : graph.getAdjacencyList()[u])
            {
                int v = edge.to;
                int weight = edge.weight;

                if (!inMST[v] && key[v] > weight)
                {
                    key[v] = weight;
                    pq.push({key[v], v});
                    parent[v] = u;
                }
            }
        }

        for (int i = 1; i < V; ++i)
        {
            mst->addEdge(parent[i], i, key[i]);
        }

        return mst;
    }
};

// Kruskal's Algorithm
class KruskalAlgorithm : public MSTAlgorithm
{
public:
    std::unique_ptr<MST> findMST(const Graph &graph) override
    {
        auto mst = std::make_unique<MST>(graph);
        int V = graph.getVerticesCount();
        std::vector<Graph::Edge> edges;

        for (int i = 0; i < V; ++i)
        {
            for (const auto &edge : graph.getAdjacencyList()[i])
            {
                if (edge.from < edge.to)
                { // To avoid duplicates in undirected graph
                    edges.push_back(edge);
                }
            }
        }

        std::sort(edges.begin(), edges.end(), [](const Graph::Edge &a, const Graph::Edge &b)
                  { return a.weight < b.weight; });

        DisjointSet ds(V);

        for (const auto &edge : edges)
        {
            int set1 = ds.find(edge.from);
            int set2 = ds.find(edge.to);

            if (set1 != set2)
            {
                mst->addEdge(edge.from, edge.to, edge.weight);
                ds.unite(set1, set2);
            }
        }

        return mst;
    }
};

// Tarjan's Algorithm (Modified for MST)
class TarjanAlgorithm : public MSTAlgorithm
{
public:
    std::unique_ptr<MST> findMST(const Graph &graph) override
    {
        int V = graph.getVerticesCount();
        std::vector<Graph::Edge> edges = gatherEdges(graph);

        // Sort all edges in ascending order
        std::sort(edges.begin(), edges.end(), [](const Graph::Edge &a, const Graph::Edge &b)
                  { return a.weight < b.weight; });

        DisjointSet ds(V);
        auto mst = std::make_unique<MST>(graph); // Pass the graph here

        for (const auto &edge : edges)
        {
            int u = ds.find(edge.from);
            int v = ds.find(edge.to);

            if (u != v)
            {
                mst->addEdge(edge.from, edge.to, edge.weight); // Use original edge vertices
                ds.unite(u, v);
            }
        }

        return mst;
    }

private:
    std::vector<Graph::Edge> gatherEdges(const Graph &graph)
    {
        std::vector<Graph::Edge> edges;
        const auto &adj = graph.getAdjacencyList();

        for (int i = 0; i < adj.size(); ++i)
        {
            for (const auto &edge : adj[i])
            {
                if (edge.from < edge.to)
                {
                    edges.push_back(edge);
                }
            }
        }

        return edges;
    }
};

class IntegerMSTAlgorithm : public MSTAlgorithm
{
};

// MST Algorithm Factory
class MSTFactory
{
public:
    enum class Algorithm
    {
        BORUVKA,
        PRIM,
        KRUSKAL,
        TARJAN,
        IntegerMST
    };

    static std::unique_ptr<MSTAlgorithm> createAlgorithm(Algorithm algo)
    {
        switch (algo)
        {
        case Algorithm::BORUVKA:
            return std::make_unique<BoruvkaAlgorithm>();
        case Algorithm::PRIM:
            return std::make_unique<PrimAlgorithm>();
        case Algorithm::KRUSKAL:
            return std::make_unique<KruskalAlgorithm>();
        case Algorithm::TARJAN:
            return std::make_unique<TarjanAlgorithm>();
        // case Algorithm::IntegerMST:
        //     return std::make_unique<IntegerMSTAlgorithm>();
        default:
            throw std::runtime_error("Unknown algorithm");
        }
    }
};

int main() {
    std::unordered_map<int, Graph> graphs;
    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream iss(line);
        std::string cmd;
        iss >> cmd;

        if (cmd == "ADD_GRAPH") {
            int graphId, numVertices;
            iss >> graphId >> numVertices;
            Graph graph(numVertices);
            int from, to, weight;
            while (iss >> from >> to >> weight) {
                graph.addEdge(from, to, weight);
            }
            graphs[graphId] = graph;
            std::cout << "Graph added successfully" << std::endl << std::endl;
        }
        else if (cmd == "SOLVE_MST") {
            int graphId;
            std::string algoName;
            iss >> graphId >> algoName;

            if (graphs.find(graphId) == graphs.end()) {
                std::cout << "Error: Graph not found" << std::endl << std::endl;
                continue;
            }

            MSTFactory::Algorithm algo;
            if (algoName == "PRIM") algo = MSTFactory::Algorithm::PRIM;
            else if (algoName == "KRUSKAL") algo = MSTFactory::Algorithm::KRUSKAL;
            else if (algoName == "BORUVKA") algo = MSTFactory::Algorithm::BORUVKA;
            else if (algoName == "TARJAN") algo = MSTFactory::Algorithm::TARJAN;
            else {
                std::cout << "Error: Unknown algorithm" << std::endl << std::endl;
                continue;
            }

            auto mstAlgorithm = MSTFactory::createAlgorithm(algo);
            auto mst = mstAlgorithm->findMST(graphs[graphId]);

            int totalWeight = mst->getTotalWeight();
            int longestDistance = mst->getLongestDistance();
            double avgDistance = mst->getAverageDistance();
            int shortestDistance = mst->getShortestEdge();

            std::cout << "Graph " << graphId << " MST Results:" << std::endl
                      << "Total Weight: " << totalWeight << std::endl
                      << "Longest Distance of MST: " << longestDistance << std::endl
                      << "Average Distance between two vertices of the graph: " << avgDistance << std::endl
                      << "Shortest Distance of MST: " << shortestDistance << "\n\n"<<std::endl; 
        }
        else {
            // std::cout << "Unknown command" << std::endl << std::endl;
        }
        std::cout << std::flush;
    }
    return 0;
}


