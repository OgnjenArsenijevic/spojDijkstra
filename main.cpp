#pragma once
#include <algorithm>
#include <assert.h>
#include <cstring>
#include <climits>
#include <functional>
#include <iostream>
#include <queue>
#include <stack>
#include <unordered_map>
#include <utility>
#include <vector>
 
template <class T>
class Graph
{
private:
	// Forward declarations
	//
	class GraphNode;
 
	// Variables
	//
	std::vector<GraphNode> vecNodes;                                // Vector to store node data.
	std::vector<std::vector<std::pair<int, long long> > > vecGraph; // Graph.
	bool m_fDirected;                                               // Whether graph is directed or not.
 
public:
	// Constructors
	//
	explicit Graph(int iArrSize, std::pair<int, int> arrEdges[], int iGraphSize, bool fDirected)
	{
		m_fDirected = fDirected;
 
		// Note that it is expected that correct values are provided for array and graph size.
		// Node ids start from 0.
		//
		for (int i = 0; i < iGraphSize; ++i)
		{
			vecGraph.push_back({});
			vecNodes.push_back(GraphNode());
		}
 
		for (int i = 0; i < iArrSize; ++i)
		{
			vecGraph[arrEdges[i].first].push_back({ arrEdges[i].second, 1LL });
 
			if (!m_fDirected)
			{
				vecGraph[arrEdges[i].second].push_back({ arrEdges[i].first, 1LL });
			}
		}
	}
 
	explicit Graph(int iArrSize, std::pair<int, int> arrEdges[], int iGraphSize, T arrData[], bool fDirected)
	{
		m_fDirected = fDirected;
 
		// Note that it is expected that correct values are provided for array and graph size.
		// Node ids start from 0.
		//
		for (int i = 0; i < iGraphSize; ++i)
		{
			vecGraph.push_back({});
			vecNodes.push_back(GraphNode(arrData[i]));
		}
 
		for (int i = 0; i < iArrSize; ++i)
		{
			vecGraph[arrEdges[i].first].push_back({ arrEdges[i].second, 1LL });
 
			if (!m_fDirected)
			{
				vecGraph[arrEdges[i].second].push_back({ arrEdges[i].first, 1LL });
			}
		}
	}
 
	explicit Graph(int iArrSize, std::pair<int, int> arrEdges[], long long arrEdgeWeights[], int iGraphSize, bool fDirected)
	{
		m_fDirected = fDirected;
 
		// Note that it is expected that correct values are provided for array and graph size.
		// Node ids start from 0.
		//
		for (int i = 0; i < iGraphSize; ++i)
		{
			vecGraph.push_back({});
			vecNodes.push_back(GraphNode());
		}
 
		for (int i = 0; i < iArrSize; ++i)
		{
			vecGraph[arrEdges[i].first].push_back({ arrEdges[i].second, arrEdgeWeights[i] });
 
			if (!m_fDirected)
			{
				vecGraph[arrEdges[i].second].push_back({ arrEdges[i].first, arrEdgeWeights[i] });
			}
		}
	}
 
	explicit Graph(int iArrSize, std::pair<int, int> arrEdges[], long long arrEdgeWeights[], int iGraphSize, T arrData[], bool fDirected)
	{
		m_fDirected = fDirected;
 
		// Note that it is expected that correct values are provided for array and graph size.
		// Node ids start from 0.
		//
		for (int i = 0; i < iGraphSize; ++i)
		{
			vecGraph.push_back({});
			vecNodes.push_back(GraphNode(arrData[i]));
		}
 
		for (int i = 0; i < iArrSize; ++i)
		{
			vecGraph[arrEdges[i].first].push_back({ arrEdges[i].second, arrEdgeWeights[i] });
 
			if (!m_fDirected)
			{
				vecGraph[arrEdges[i].second].push_back({ arrEdges[i].first, arrEdgeWeights[i] });
			}
		}
	}
 
	explicit Graph(std::vector<std::pair<int, int> > vecEdges, int iGraphSize, bool fDirected)
	{
		m_fDirected = fDirected;
 
		// Note that it is expected that correct values are provided for array and graph size.
		// Node ids start from 0.
		//
		for (int i = 0; i < iGraphSize; ++i)
		{
			vecGraph.push_back({});
			vecNodes.push_back(GraphNode());
		}
 
		for (unsigned int i = 0; i < vecEdges.size(); ++i)
		{
			vecGraph[vecEdges[i].first].push_back({ vecEdges[i].second, 1LL });
 
			if (!m_fDirected)
			{
				vecGraph[vecEdges[i].second].push_back({ vecEdges[i].first, 1LL });
			}
		}
	}
 
	explicit Graph(std::vector<std::pair<int, int> > vecEdges, int iGraphSize, std::vector<T> vecData, bool fDirected)
	{
		m_fDirected = fDirected;
 
		// Note that it is expected that correct values are provided for array and graph size.
		// Node ids start from 0.
		//
		for (int i = 0; i < iGraphSize; ++i)
		{
			vecGraph.push_back({});
			vecNodes.push_back(GraphNode(vecData[i]));
		}
 
		for (unsigned int i = 0; i < vecEdges.size(); ++i)
		{
			vecGraph[vecEdges[i].first].push_back({ vecEdges[i].second, 1LL });
 
			if (!m_fDirected)
			{
				vecGraph[vecEdges[i].second].push_back({ vecEdges[i].first, 1LL });
			}
		}
	}
 
	explicit Graph(std::vector<std::pair<int, int> > vecEdges, std::vector<long long> vecEdgeWeights, int iGraphSize, bool fDirected)
	{
		m_fDirected = fDirected;
 
		// Note that it is expected that correct values are provided for array and graph size.
		// Node ids start from 0.
		//
		for (int i = 0; i < iGraphSize; ++i)
		{
			vecGraph.push_back({});
			vecNodes.push_back(GraphNode());
		}
 
		for (unsigned int i = 0; i < vecEdges.size(); ++i)
		{
			vecGraph[vecEdges[i].first].push_back({ vecEdges[i].second, vecEdgeWeights[i] });
 
			if (!m_fDirected)
			{
				vecGraph[vecEdges[i].second].push_back({ vecEdges[i].first, vecEdgeWeights[i] });
			}
		}
	}
 
	explicit Graph(std::vector<std::pair<int, int> > vecEdges, std::vector<long long> vecEdgeWeights, int iGraphSize, std::vector<T> vecData, bool fDirected)
	{
		m_fDirected = fDirected;
 
		// Note that it is expected that correct values are provided for array and graph size.
		// Node ids start from 0.
		//
		for (int i = 0; i < iGraphSize; ++i)
		{
			vecGraph.push_back({});
			vecNodes.push_back(GraphNode(vecData[i]));
		}
 
		for (unsigned int i = 0; i < vecEdges.size(); ++i)
		{
			vecGraph[vecEdges[i].first].push_back({ vecEdges[i].second, vecEdgeWeights[i] });
 
			if (!m_fDirected)
			{
				vecGraph[vecEdges[i].second].push_back({ vecEdges[i].first, vecEdgeWeights[i] });
			}
		}
	}
 
	// Methods
	//
	std::vector<bool> depth_first_search(const int& iNodeId, std::vector<bool> vecVisited = std::vector<bool>());					// Does depth first search starting from node with given id.
	std::vector<long long> breadth_first_search(const int& iNodeId, std::vector<long long> vecDistance = std::vector<long long>());	// Does breadth first search starting from node with given id.
	std::vector<long long> dijkstra(const int& iNodeId, std::vector<long long> vecDistance = std::vector<long long>());				// Does Dijkstra algorithm starting from node with given id.
	std::vector<long long> bellman_ford(const int& iNodeId, std::vector<long long> vecDistance = std::vector<long long>());			// Does Bellman-Ford algorithm starting from node with given id.
	std::vector<std::vector<long long> > floyd_warshall();																			// Does Floyd-Warshall alogrithm to find shortest distances between every pair of nodes.
	std::pair<long long, std::vector<std::pair<long long, std::pair<int, int> > > > prim();											// Does Prim algorithm to find minimum spanning tree of the graph.
	std::pair<long long, std::vector<std::pair<long long, std::pair<int, int> > > > kruskal();										// Does Kruskal algorithm to find minimum spanning tree of the graph.
	std::pair<int, std::vector<std::vector<int> > > connected_compontents();														// Returns number of connected components in the graph and ids of the nodes in each component.
	std::pair<int, std::vector<std::vector<int> > > kosaraju();																		// Returns number of strongly connected components in the graph and ids of the nodes in each component.
	std::pair<int, std::vector<std::vector<int> > > tarjan();																		// Returns number of strongly connected components in the graph and ids of the nodes in each component.
	std::pair<int, std::vector<int> > articulation_points();																		// Returns articulation points of the graph.
	std::pair<int, std::vector<std::pair<int, int>> > bridges();																	// Returns bridges of the graph.
	long long edmonds_karp(const int& iSourceNodeId, const int& iSinkNodeId);														// Returns max flow from source to sink using Edmonds-Karp algorithm.
	long long dinic(const int& iSourceNodeId, const int& iSinkNodeId);																// Returns max flow from source to sink using Dinics algorithm.
	bool cyclic();																													// Returns whether graph has cycle.
 
private:
	class GraphNode
	{
	private:
		friend class Graph;		// Graph class methods need to access GraphNode info.
		T m_xData;				// Node data.
 
		// Constructors
		//
		GraphNode() : m_xData(0) {};
 
		GraphNode(T xData) : m_xData(xData) {};
	};
};
 
//-------------------------------------------------------------------------------------------
// Function: depth_first_search
//
// Description:
//	Does depth first search from node with given id.
//
// Parameters:
//	iNodeId - Id of the node from which dfs should be started.
//
// Note:
//	It expected that node id is valid. In case that node id is
//  not valid this function will assert.
//
// Time complexity:
//	O(V + E) - where V is the number of vertices and E is the number of edges in the graph.
//
template<class T>
std::vector<bool> Graph<T>::depth_first_search(const int& iNodeId, std::vector<bool> vecVisited)
{
	if ((unsigned int)iNodeId >= vecGraph.size())
	{
		// Invalid node id.
		//
		assert(false);
	}
 
	std::vector<bool> vecResultVisited;
 
	if (vecVisited.size() == vecGraph.size())
	{
		vecResultVisited = vecVisited;
	}
	else
	{
		vecResultVisited.insert(vecResultVisited.begin(), vecGraph.size(), false);
	}
 
	auto dfs = [&](auto&& dfs, const int& iNodeId) -> void
	{
		vecResultVisited[iNodeId] = true;
 
		for (auto& it : vecGraph[iNodeId])
		{
			int iNeighbourNodeId = it.first;
 
			if (!vecResultVisited[iNeighbourNodeId])
			{
				dfs(dfs, iNeighbourNodeId);
			}
		}
	};
 
	dfs(dfs, iNodeId);
 
	return vecResultVisited;
}
 
//-------------------------------------------------------------------------------------------
// Function: breadth_first_search
//
// Description:
//	Does breadth first search from node with given id.
//
// Parameters:
//	iNodeId - Id of the node from which bfs should be started.
//
// Note:
//	It expected that node id is valid. In case that node id is
//  not valid this function will assert.
//
// Time complexity:
//	O(V + E) - where V is the number of vertices and E is the number of edges in the graph.
//
template<class T>
std::vector<long long> Graph<T>::breadth_first_search(const int& iNodeId, std::vector<long long> vecDistance)
{
	std::vector<long long> vecResultDistance;
 
	if (vecDistance.size() == vecGraph.size())
	{
		vecResultDistance = vecDistance;
	}
	else
	{
		vecResultDistance.insert(vecResultDistance.begin(), vecGraph.size(), 0LL);
	}
 
	std::queue<int> que;
	que.push(iNodeId);
 
	while (!que.empty())
	{
		int iCurrentNodeId = que.front();
		que.pop();
 
		for (auto& it : vecGraph[iCurrentNodeId])
		{
			int iNeighbourNodeId = it.first;
 
			if (!vecResultDistance[iNeighbourNodeId] != 0)
			{
				vecResultDistance[iNeighbourNodeId] = vecResultDistance[iCurrentNodeId] + 1LL;
				que.push(iNeighbourNodeId);
			}
		}
	}
 
	return vecResultDistance;
}
 
//-------------------------------------------------------------------------------------------
// Function: dijkstra
//
// Description:
//	Does Dijkstra algorithm from node with given id.
//
// Parameters:
//	iNodeId - Id of the node from which bfs should be started.
//
// Note:
//	It expected that node id is valid. In case that node id is
//  not valid this function will assert.
//  Is is also expected that all edges have positive weight,
//  otherwise function will assert.
//
// Time complexity:
//	O(ElogV) - where V is the number of vertices and E is the number of edges in the graph.
//
template<class T>
std::vector<long long> Graph<T>::dijkstra(const int& iNodeId, std::vector<long long> vecDistance)
{
	if ((unsigned int)iNodeId >= vecGraph.size())
	{
		// Invalid node id.
		//
		assert(false);
	}
 
	std::vector<long long> vecResultDistance;
 
	if (vecDistance.size() == vecGraph.size())
	{
		vecResultDistance = vecDistance;
	}
	else
	{
		vecResultDistance.insert(vecResultDistance.begin(), vecGraph.size(), LLONG_MAX);
	}
 
	std::priority_queue<std::pair<long long, int>, std::vector<std::pair<long long,  int> >, std::greater<std::pair<long long, int> > > pq;
	vecResultDistance[iNodeId] = 0;
	pq.push({ 0, iNodeId });
 
	while (!pq.empty())
	{
		long long llCurrentDistance = pq.top().first;
		int iCurrentNodeId = pq.top().second;
		pq.pop();
 
		if (llCurrentDistance > vecResultDistance[iCurrentNodeId])
		{
			continue;
		}
 
		for (auto& it : vecGraph[iCurrentNodeId])
		{
			int iNeighbourNodeId = it.first;
			long long llWeight = it.second;
 
			if (llWeight <= 0LL)
			{
				// All edges must have positive weight.
				//
				assert(false);
			}
 
			if (vecResultDistance[iCurrentNodeId] + llWeight < vecResultDistance[iNeighbourNodeId])
			{
				vecResultDistance[iNeighbourNodeId] = vecResultDistance[iCurrentNodeId] + llWeight;
				pq.push({ vecResultDistance[iNeighbourNodeId], iNeighbourNodeId });
			}
		}
	}
 
	return vecResultDistance;
}
 
//-------------------------------------------------------------------------------------------
// Function: bellman_ford
//
// Description:
//	Does Bellman-Ford algorithm from node with given id.
//
// Parameters:
//	iNodeId - Id of the node from which bfs should be started.
//
// Note:
//	It expected that node id is valid. In case that node id is
//  not valid this function will assert.
//  In case that graph contains negative weight cycle function
//  will assert.
//
// Time complexity:
//	O(VE) - where V is the number of vertices and E is the number of edges in the graph.
//
template<class T>
std::vector<long long> Graph<T>::bellman_ford(const int& iNodeId, std::vector<long long> vecDistance)
{
	if ((unsigned int)iNodeId >= vecGraph.size())
	{
		// Invalid node id.
		//
		assert(false);
	}
 
	std::vector<long long> vecResultDistance;
 
	if (vecDistance.size() == vecGraph.size())
	{
		vecResultDistance = vecDistance;
	}
	else
	{
		vecResultDistance.insert(vecResultDistance.begin(), vecGraph.size(), LLONG_MAX);
	}
 
	vecResultDistance[iNodeId] = 0;
 
	for (unsigned int i = 0; i < vecGraph.size() - 1; ++i)
	{
		for (unsigned int j = 0; j < vecGraph.size(); ++j)
		{
			for (auto& it : vecGraph[j])
			{
				int iNeighbourNodeId = it.first;
				long long llEdgeWeight = it.second;
				vecResultDistance[iNeighbourNodeId] = std::min(vecResultDistance[iNeighbourNodeId], vecResultDistance[j] + llEdgeWeight);
			}
		}
	}
 
	for (unsigned int i = 0; i < vecGraph.size(); ++i)
	{
		for (auto& it : vecGraph[i])
		{
			int iNeighbourNodeId = it.first;
			long long llEdgeWeight = it.second;
 
			if (vecResultDistance[i] != LLONG_MAX && vecResultDistance[i] + llEdgeWeight < vecResultDistance[iNeighbourNodeId])
			{
				// Graph contains negative edge cycle.
				//
				assert(false);
			}
		}
	}
 
	return vecResultDistance;
}
 
//-------------------------------------------------------------------------------------------
// Function: floyd_warshall
//
// Description:
//  Does Floyd-Warshall alogrithm to find shortest distances between every pair of nodes.
//
// Time complexity:
//	O(V^3) - where V is the number of vertices and E is the number of edges in the graph.
//
template<class T>
std::vector<std::vector<long long> > Graph<T>::floyd_warshall()
{
	std::vector<std::vector<long long> > vecResult(vecGraph.size(), std::vector<long long>(vecGraph.size(), LLONG_MAX));
 
	for (unsigned int i = 0; i < vecGraph.size(); ++i)
	{
		for (auto& it : vecGraph[i])
		{
			int iNeighbourNodeId = it.first;
			long long llEdgeWeight = it.second;
			vecResult[i][iNeighbourNodeId] = llEdgeWeight;
		}
 
		vecResult[i][i] = 0LL;
	}
 
	for (unsigned int k = 0; k < vecGraph.size(); ++k)
	{
		for (unsigned int i = 0; i < vecGraph.size(); ++i)
		{
			for (unsigned int j = 0; j < vecGraph.size(); ++j)
			{
				if (vecResult[i][j] > vecResult[i][k] + vecResult[k][j] && vecResult[k][j] != LLONG_MAX && vecResult[i][k] != LLONG_MAX)
				{
					vecResult[i][j] = vecResult[i][k] + vecResult[k][j];
				}
			}
		}
	}
 
	return vecResult;
}
 
//-------------------------------------------------------------------------------------------
// Function: prim
//
// Description:
//	Does Prim algorithm to find minumum spanning tree in the graph.
//
// Note:
//	This function works only for undirected graph. In case that graph is directed
//  function will assert.
//
// Time complexity:
//	O(ElogV) - where V is the number of vertices and E is the number of edges in the graph.
//
template<class T>
std::pair<long long, std::vector<std::pair<long long, std::pair<int, int> > > > Graph<T>::prim()
{
	if (m_fDirected)
	{
		// Graph must be undirected.
		//
		assert(false);
	}
 
	std::vector<bool> vecUsed(vecGraph.size(), false);
	std::priority_queue<std::pair<long long, std::pair<int, int> >, std::vector<std::pair<long long, std::pair<int, int> > >, std::greater<std::pair<long long, std::pair<int, int> > > > pq;
 
	auto process = [&](const int& iNodeId) -> void
	{
		vecUsed[iNodeId] = true;
 
		for (auto& it : vecGraph[iNodeId])
		{
			int iNeighbourNodeId = it.first;
			long long llEdgeWeight = it.second;
 
			if (!vecUsed[iNeighbourNodeId])
			{
				pq.push({ llEdgeWeight, {iNodeId, iNeighbourNodeId} });
			}
		}
	};
 
	long long llResult = 0LL;
	std::vector<std::pair<long long, std::pair<int, int> > > vecResult;
 
	process(0);
	while (!pq.empty())
	{
		int iNodeId1 = pq.top().second.first;
		int iNodeId2 = pq.top().second.second;
		long long llEdgeWeight = pq.top().first;
		pq.pop();
 
		if (!vecUsed[iNodeId2])
		{
			vecResult.push_back({ llEdgeWeight, {iNodeId1, iNodeId2} });
			llResult += llEdgeWeight;
			process(iNodeId2);
		}
	}
 
	return { llResult, vecResult };
}
 
//-------------------------------------------------------------------------------------------
// Function: kruskal
//
// Description:
//	Does Kruskal algorithm to find minumum spanning tree in the graph.
//
// Note:
//	This function works only for undirected graph. In case that graph is directed
//  function will assert.
//
// Time complexity:
//	O(ElogV) - where V is the number of vertices and E is the number of edges in the graph.
//
template<class T>
std::pair<long long, std::vector<std::pair<long long, std::pair<int, int> > > > Graph<T>::kruskal()
{
	if (m_fDirected)
	{
		// Graph must be undirected.
		//
		assert(false);
	}
 
	std::vector<int> vecParent(vecGraph.size());
	std::vector<int> vecRank(vecGraph.size(), 0);
 
	for (unsigned int i = 0; i < vecGraph.size(); ++i)
	{
		vecParent[i] = i;
	}
 
	auto find = [&](auto&& find, const int& iNodeId) -> int
	{
		if (vecParent[iNodeId] == iNodeId)
		{
			return vecParent[iNodeId];
		}
 
		vecParent[iNodeId] = find(find, vecParent[iNodeId]);
		return vecParent[iNodeId];
	};
 
	auto is_in_same_set = [&](const int& iNodeId1, const int& iNodeId2) -> bool
	{
		return find(find, iNodeId1) == find(find, iNodeId2);
	};
 
	auto unite = [&](const int& iNodeId1, const int& iNodeId2) -> void
	{
		int iParentNodeId1 = find(find, iNodeId1);
		int iParentNodeId2 = find(find, iNodeId2);
 
		if (iParentNodeId1 == iParentNodeId2)
		{
			return;
		}
 
		if (vecRank[iParentNodeId1] < vecRank[iParentNodeId2])
		{
			vecParent[iParentNodeId1] = iParentNodeId2;
		}
		else
		{
			vecParent[iParentNodeId2] = iParentNodeId1;
 
			if (vecRank[iParentNodeId1] == vecRank[iParentNodeId2])
			{
				++vecRank[iParentNodeId1];
			}
		}
	};
 
	std::vector<std::pair<long long, std::pair<int, int> > > vecEdges;
 
	for (unsigned int i = 0; i < vecGraph.size(); ++i)
	{
		for (auto& it : vecGraph[i])
		{
			int iNeighbourNodeId = it.first;
			long long llEdgeWeight = it.second;
			vecEdges.push_back({ llEdgeWeight, {i, iNeighbourNodeId} });
		}
	}
 
	std::sort(vecEdges.begin(), vecEdges.end());
 
	long long llResult = 0LL;
	std::vector<std::pair<long long, std::pair<int, int> > > vecResult;
 
	for (unsigned int i = 0; i < vecEdges.size(); ++i)
	{
		int iNodeId1 = vecEdges[i].second.first;
		int iNodeId2 = vecEdges[i].second.second;
		long long llEdgeWeight = vecEdges[i].first;
 
		if (!is_in_same_set(iNodeId1, iNodeId2))
		{
			vecResult.push_back({ llEdgeWeight, {iNodeId1, iNodeId2} });
			llResult += llEdgeWeight;
			unite(iNodeId1, iNodeId2);
		}
	}
 
	return { llResult, vecResult };
}
 
//-------------------------------------------------------------------------------------------
// Function: connected_compontents
//
// Description:
//  Returns number of connected components in the graph and ids of the nodes in each component.
//
// Note:
//	It expected that graph is undirected, otherwise function will assert.
//
// Time complexity:
//	O(V + E) - where V is the number of vertices and E is the number of edges in the graph.
//
template<class T>
 std::pair<int, std::vector<std::vector<int> > > Graph<T>::connected_compontents()
{
	 if (m_fDirected)
	 {
		 // Graph must be undirected.
		 //
		 assert(false);
	 }
 
	 std::vector<bool> vecVisited(vecGraph.size(), false);
 
	 auto dfs = [&](auto&& dfs, const int& iNodeId, std::vector<int>& vecCurrentComponent) -> void
	 {
		 vecCurrentComponent.push_back(iNodeId);
		 vecVisited[iNodeId] = true;
 
		 for (auto& it : vecGraph[iNodeId])
		 {
			 int iNeighbourNodeId = it.first;
 
			 if (!vecVisited[iNeighbourNodeId])
			 {
				 dfs(dfs, iNeighbourNodeId, vecCurrentComponent);
			 }
		 }
	 };
 
	 int iResult = 0;
	 std::vector<std::vector<int> > vecResult;
 
	 for (unsigned int i = 0; i < vecGraph.size(); ++i)
	 {
		 if (!vecVisited[i])
		 {
			 std::vector<int> vecCurrentComponent;
			 dfs(dfs, i, vecCurrentComponent);
			 ++iResult;
			 vecResult.push_back(vecCurrentComponent);
		 }
	 }
 
	 return { iResult, vecResult };
}
 
//-------------------------------------------------------------------------------------------
// Function: kosaraju
//
// Description:
//	Returns number of strongly connected components in the graph and ids of the nodes in each component.
//
// Note:
//	It expected that graph is directed, otherwise function will assert.
//
// Time complexity:
//	O(V + E) - where V is the number of vertices and E is the number of edges in the graph.
//
template<class T>
std::pair<int, std::vector<std::vector<int> > > Graph<T>::kosaraju()
{
	if (!m_fDirected)
	{
		// Graph must be directed.
		//
		assert(false);
	}
 
	std::stack<int> st;
	std::vector<bool> vecVisited(vecGraph.size(), false);
 
	auto dfs1 = [&](auto&& dfs1, const int& iNodeId) -> void
	{
		vecVisited[iNodeId] = true;
 
		for (auto& it : vecGraph[iNodeId])
		{
			int iNeighbourNodeId = it.first;
 
			if (!vecVisited[iNeighbourNodeId])
			{
				dfs1(dfs1, iNeighbourNodeId);
			}
		}
 
		st.push(iNodeId);
	};
 
	std::vector<std::vector<int>> vecTransponseGraph(vecGraph.size());
 
	auto dfs2 = [&](auto&& dfs2, const int& iNodeId, std::vector<int>& vecCurrentComponent) -> void
	{
		vecCurrentComponent.push_back(iNodeId);
		vecVisited[iNodeId] = true;
 
		for (auto& it : vecTransponseGraph[iNodeId])
		{
			int iNeighbourNodeId = it;
 
			if (!vecVisited[iNeighbourNodeId])
			{
				dfs2(dfs2, iNeighbourNodeId, vecCurrentComponent);
			}
		}
	};
 
	for (unsigned int i = 0; i < vecGraph.size(); ++i)
	{
		if (!vecVisited[i])
		{
			dfs1(dfs1, i);
		}
	}
 
	// Make transponse graph.
	//
	for (unsigned int i = 0; i < vecGraph.size(); ++i)
	{
		for (auto& it : vecGraph[i])
		{
			vecTransponseGraph[it.first].push_back(i);
		}
 
		vecVisited[i] = false;
	}
 
	int iResult = 0;
	std::vector<std::vector<int> > vecResult;
 
	while (!st.empty())
	{
		int iCurrentNode = st.top();
		st.pop();
 
		if (!vecVisited[iCurrentNode])
		{
			std::vector<int> vecCurrentComponent;
			dfs2(dfs2, iCurrentNode, vecCurrentComponent);
			++iResult;
			vecResult.push_back(vecCurrentComponent);
		}
	}
 
	return { iResult, vecResult };
}
 
//-------------------------------------------------------------------------------------------
// Function: tarjan
//
// Description:
//	Returns number of strongly connected components in the graph and ids of the nodes in each component.
//
// Note:
//	It expected that graph is directed, otherwise function will assert.
//
// Time complexity:
//	O(V + E) - where V is the number of vertices and E is the number of edges in the graph.
//
template<class T>
std::pair<int, std::vector<std::vector<int> > > Graph<T>::tarjan()
{
	if (!m_fDirected)
	{
		// Graph must be directed.
		//
		assert(false);
	}
 
	int iTime = 0;
	std::vector<int> vecDiscovery(vecGraph.size(), -1);
	std::vector<int> vecLow(vecGraph.size(), -1);
	std::vector<bool> vecOnStack(vecGraph.size(), false);
	std::stack<int> st;
 
	int iResult = 0;
	std::vector<std::vector<int> > vecResult;
 
	auto dfs = [&](auto&& dfs, const int& iNodeId) -> void
	{
		vecDiscovery[iNodeId] = iTime;
		vecLow[iNodeId] = iTime;
		++iTime;
		st.push(iNodeId);
		vecOnStack[iNodeId] = true;
 
		for (auto& it : vecGraph[iNodeId])
		{
			int iNeighbourNodeId = it.first;
 
			if (vecDiscovery[iNeighbourNodeId] == -1)
			{
				dfs(dfs, iNeighbourNodeId);
 
				vecLow[iNodeId] = std::min(vecLow[iNodeId], vecLow[iNeighbourNodeId]);
			}
			else if (vecOnStack[iNeighbourNodeId])
			{
				vecLow[iNodeId] = std::min(vecLow[iNodeId], vecDiscovery[iNeighbourNodeId]);
			}
		}
 
		if (vecLow[iNodeId] == vecDiscovery[iNodeId])
		{
			std::vector<int> vecCurrentComponent;
 
			while (st.top() != iNodeId)
			{
				vecCurrentComponent.push_back(st.top());
				vecOnStack[st.top()] = false;
				st.pop();
			}
 
			vecCurrentComponent.push_back(st.top());
			vecOnStack[st.top()] = false;
			st.pop();
			vecResult.push_back(vecCurrentComponent);
			++iResult;
		}
	};
 
	for (unsigned int i = 0; i < vecGraph.size(); ++i)
	{
		if (vecDiscovery[i] == -1)
		{
			dfs(dfs, i);
		}
	}
 
	return { iResult, vecResult };
}
 
//-------------------------------------------------------------------------------------------
// Function: articulation_points
//
// Description:
//	Returns articulation points of the graph.
//
// Time complexity:
//	O(V + E) - where V is the number of vertices and E is the number of edges in the graph.
//
template<class T>
std::pair<int, std::vector<int> > Graph<T>::articulation_points()
{
	int iTime = 0;
	std::vector<int> vecDiscovery(vecGraph.size(), -1);
	std::vector<int> vecLow(vecGraph.size(), -1);
	std::vector<bool> vecVisited(vecGraph.size(), false);
 
	int iResult = 0;
	std::vector<int> vecResult;
 
	auto dfs = [&](auto&& dfs, const int& iNodeId, const int& iParentNodeId) -> void
	{
		int iChildren = 0;
		vecVisited[iNodeId] = true;
		vecDiscovery[iNodeId] = iTime;
		vecLow[iNodeId] = iTime;
		++iTime;
 
		for (auto& it : vecGraph[iNodeId])
		{
			int iNeighbourNodeId = it.first;
 
			if (!vecVisited[iNeighbourNodeId])
			{
				++iChildren;
				dfs(dfs, iNeighbourNodeId, iNodeId);
 
				vecLow[iNodeId] = std::min(vecLow[iNodeId], vecLow[iNeighbourNodeId]);
 
				if (iParentNodeId != -1 && vecLow[iNeighbourNodeId] >= vecDiscovery[iNodeId])
				{
					vecResult.push_back(iNodeId);
					++iResult;
				}
			}
			else if (iNeighbourNodeId != iParentNodeId)
			{
				vecLow[iNodeId] = std::min(vecLow[iNodeId], vecDiscovery[iNeighbourNodeId]);
			}
		}
 
		if (iParentNodeId == -1 && iChildren > 1)
		{
			vecResult.push_back(iNodeId);
			++iResult;
		}
	};
 
	for (unsigned int i = 0; i < vecGraph.size(); ++i)
	{
		if (!vecVisited[i])
		{
			dfs(dfs, i, -1);
		}
	}
 
	return { iResult, vecResult };
}
 
//-------------------------------------------------------------------------------------------
// Function: bridges
//
// Description:
//	Returns bridges of the graph.
//
// Time complexity:
//	O(V + E) - where V is the number of vertices and E is the number of edges in the graph.
//
template<class T>
std::pair<int, std::vector<std::pair<int, int> > > Graph<T>::bridges()
{
	int iTime = 0;
	std::vector<int> vecDiscovery(vecGraph.size(), -1);
	std::vector<int> vecLow(vecGraph.size(), -1);
	std::vector<bool> vecVisited(vecGraph.size(), false);
 
	int iResult = 0;
	std::vector<std::pair<int, int> > vecResult;
 
	auto dfs = [&](auto&& dfs, const int& iNodeId, const int& iParentNodeId) -> void
	{
		vecVisited[iNodeId] = true;
		vecDiscovery[iNodeId] = iTime;
		vecLow[iNodeId] = iTime;
		++iTime;
 
		for (auto& it : vecGraph[iNodeId])
		{
			int iNeighbourNodeId = it.first;
 
			if (!vecVisited[iNeighbourNodeId])
			{
				dfs(dfs, iNeighbourNodeId, iNodeId);
 
				vecLow[iNodeId] = std::min(vecLow[iNodeId], vecLow[iNeighbourNodeId]);
 
				if (vecLow[iNeighbourNodeId] > vecDiscovery[iNodeId])
				{
					vecResult.push_back({iNodeId, iNeighbourNodeId});
					++iResult;
				}
			}
			else if (iNeighbourNodeId != iParentNodeId)
			{
				vecLow[iNodeId] = std::min(vecLow[iNodeId], vecDiscovery[iNeighbourNodeId]);
			}
		}
	};
 
	for (unsigned int i = 0; i < vecGraph.size(); ++i)
	{
		if (!vecVisited[i])
		{
			dfs(dfs, i, -1);
		}
	}
 
	return { iResult, vecResult };
}
 
//-------------------------------------------------------------------------------------------
// Function: edmonds_karp
//
// Description:
//	Returns max flow from source to sink using Edmonds-Karp algorithm.
//
// Parameters:
//	iSourceNodeId - Id of the source node.
//	iSinkNodeId - Id of the sink node.
//
// Note:
//	It expected that node ids are valid. In case that node id is
//  not valid this function will assert.
//
// Time complexity:
//	O(VE^2) - where V is the number of vertices and E is the number of edges in the graph.
//
template<class T>
long long Graph<T>::edmonds_karp(const int& iSourceNodeId, const int& iSinkNodeId)
{
	std::vector<std::vector<long long> > vecCapacity(vecGraph.size(), std::vector<long long>(vecGraph.size(), 0LL));
	std::vector<std::vector<long long> > vecPassedFlow(vecGraph.size(), std::vector<long long>(vecGraph.size(), 0LL));
	long long llResult = 0LL;
 
	std::vector<std::vector<int> > vecResidualGraph;
 
	for (unsigned int i = 0; i < vecGraph.size(); ++i)
	{
		std::vector<int> vecCurrent;
 
		for (auto& it : vecGraph[i])
		{
			vecCurrent.push_back(it.first);
		}
 
		vecResidualGraph.push_back(vecCurrent);
	}
 
	if (m_fDirected)
	{
		for (unsigned int i = 0; i < vecGraph.size(); ++i)
		{
 
			for (auto& it : vecGraph[i])
			{
				vecResidualGraph[it.first].push_back(i);
			}
		}
	}
 
 
	for (unsigned int i = 0; i < vecGraph.size(); ++i)
	{
		for (auto& it : vecGraph[i])
		{
			int iNeighbourNodeId = it.first;
			long long llEdgeWeight = it.second;
			vecCapacity[i][iNeighbourNodeId] = llEdgeWeight;
		}
	}
 
	while (true)
	{
		long long llCurrentFlow = 0;
		std::vector<int> vecParent(vecGraph.size(), -1);
		std::vector<long long> vecCurrPathCapacity(vecGraph.size(), 0LL);
		std::queue<int> que;
		vecParent[iSourceNodeId] = -2;
		vecCurrPathCapacity[iSourceNodeId] = LLONG_MAX;
		que.push(iSourceNodeId);
 
		while (!que.empty())
		{
			int iCurrentNodeId = que.front();
			que.pop();
 
			for (auto& it : vecResidualGraph[iCurrentNodeId])
			{
				int iNeighbourNodeId = it;
				if (vecParent[iNeighbourNodeId] == -1 && vecCapacity[iCurrentNodeId][iNeighbourNodeId] - vecPassedFlow[iCurrentNodeId][iNeighbourNodeId] > 0)
				{
					vecParent[iNeighbourNodeId] = iCurrentNodeId;
					vecCurrPathCapacity[iNeighbourNodeId] = std::min(vecCurrPathCapacity[iCurrentNodeId], vecCapacity[iCurrentNodeId][iNeighbourNodeId] - vecPassedFlow[iCurrentNodeId][iNeighbourNodeId]);
 
					if (iNeighbourNodeId == iSinkNodeId)
					{
						llCurrentFlow = vecCurrPathCapacity[iNeighbourNodeId];
						break;
					}
					que.push(iNeighbourNodeId);
				}
			}
		}
 
		if (llCurrentFlow == 0LL)
		{
			break;
		}
 
		llResult += llCurrentFlow;
 
		int iCurrentNodeId = iSinkNodeId;
 
		while (iCurrentNodeId != iSourceNodeId)
		{
			vecCapacity[vecParent[iCurrentNodeId]][iCurrentNodeId] -= llCurrentFlow;
			vecCapacity[iCurrentNodeId][vecParent[iCurrentNodeId]] += llCurrentFlow;
			iCurrentNodeId = vecParent[iCurrentNodeId];
		}
	}
 
	return llResult;
}
 
//-------------------------------------------------------------------------------------------
// Function: dinic
//
// Description:
//	Returns max flow from source to sink using Edmonds-Karp algorithm.
//
// Parameters:
//	iSourceNodeId - Id of the source node.
//	iSinkNodeId - Id of the sink node.
//
// Note:
//	It expected that node ids are valid. In case that node id is
//  not valid this function will assert.
//
// Time complexity:
//	O(EV^2) - where V is the number of vertices and E is the number of edges in the graph.
//
template<class T>
long long Graph<T>::dinic(const int& iSourceNodeId, const int& iSinkNodeId)
{
	std::vector<std::vector<long long> > vecCapacity(vecGraph.size(), std::vector<long long>(vecGraph.size(), 0LL));
	std::vector<std::vector<long long> > vecPassedFlow(vecGraph.size(), std::vector<long long>(vecGraph.size(), 0LL));
	std::vector<int> vecLevel(vecGraph.size());
	long long llResult = 0LL;
 
	std::vector<std::vector<int> > vecResidualGraph;
 
	for (unsigned int i = 0; i < vecGraph.size(); ++i)
	{
		std::vector<int> vecCurrent;
 
		for (auto& it : vecGraph[i])
		{
			vecCurrent.push_back(it.first);
		}
 
		vecResidualGraph.push_back(vecCurrent);
	}
 
	if (m_fDirected)
	{
		for (unsigned int i = 0; i < vecGraph.size(); ++i)
		{
 
			for (auto& it : vecGraph[i])
			{
				vecResidualGraph[it.first].push_back(i);
			}
		}
	}
 
 
	for (unsigned int i = 0; i < vecGraph.size(); ++i)
	{
		for (auto& it : vecGraph[i])
		{
			int iNeighbourNodeId = it.first;
			long long llEdgeWeight = it.second;
			vecCapacity[i][iNeighbourNodeId] = llEdgeWeight;
		}
	}
 
	auto bfs = [&]() -> bool
	{
		for (unsigned int i = 0; i < vecGraph.size(); ++i)
		{
			vecLevel[i] = -1;
		}
 
		vecLevel[iSourceNodeId] = 0;
		std::queue<int> que;
		que.push(iSourceNodeId);
 
		while (!que.empty())
		{
			int iCurrentNodeId = que.front();
			que.pop();
 
			for (auto& it : vecResidualGraph[iCurrentNodeId])
			{
				int iNeighbourNodeId = it;
				if (vecLevel[iNeighbourNodeId] < 0 && vecCapacity[iCurrentNodeId][iNeighbourNodeId] - vecPassedFlow[iCurrentNodeId][iNeighbourNodeId] > 0)
				{
					vecLevel[iNeighbourNodeId] = vecLevel[iCurrentNodeId] + 1;
					que.push(iNeighbourNodeId);
				}
			}
		}
 
		return vecLevel[iSinkNodeId] < 0 ? false : true;
	};
 
	auto dfs = [&](auto&& dfs, const int& iCurrentNodeId, long long llFlow, std::vector<int>& vecStart) -> long long
	{
		if (iCurrentNodeId == iSinkNodeId)
		{
			return llFlow;
		}
 
		while ((unsigned int)vecStart[iCurrentNodeId] < vecResidualGraph[iCurrentNodeId].size())
		{
			int iNeighbourNodeId = vecResidualGraph[iCurrentNodeId][vecStart[iCurrentNodeId]];
 
			if (vecLevel[iNeighbourNodeId] == vecLevel[iCurrentNodeId] + 1 && vecCapacity[iCurrentNodeId][iNeighbourNodeId] - vecPassedFlow[iCurrentNodeId][iNeighbourNodeId] > 0)
			{
				long long llCurrFlow = std::min(llFlow, vecCapacity[iCurrentNodeId][iNeighbourNodeId] - vecPassedFlow[iCurrentNodeId][iNeighbourNodeId]);
				long long llTempFlow = dfs(dfs, iNeighbourNodeId, llCurrFlow, vecStart);
 
				if (llTempFlow > 0)
				{
					vecPassedFlow[iCurrentNodeId][iNeighbourNodeId] += llTempFlow;
					vecPassedFlow[iNeighbourNodeId][iCurrentNodeId] -= llTempFlow;
					return llTempFlow;
				}
			}
 
			++vecStart[iCurrentNodeId];
		}
 
		return 0LL;
	};
 
	if (iSourceNodeId == iSinkNodeId)
	{
		return -1LL;
	}
 
	while (bfs())
	{
		std::vector<int> vecStart(vecGraph.size(), 0);
 
		while (long long llFlow = dfs(dfs, iSourceNodeId, LLONG_MAX, vecStart))
		{
			llResult += llFlow;
		}
	}
 
	return llResult;
}
 
//-------------------------------------------------------------------------------------------
// Function: cyclic
//
// Description:
//	Returns whether graph has cycle.
//
// Parameters:
//
// Time complexity:
//	O(V + E) - where V is the number of vertices and E is the number of edges in the graph.
//
template<class T>
bool Graph<T>::cyclic()
{
	const int WHITE = 0;
	const int GRAY = 1;
	const int BLACK = 2;
	std::vector<int> vecVisited (vecGraph.size(), 0);
 
	auto cyclic_directed = [&](auto&& cyclic_directed, int iNodeId) -> bool
	{
		vecVisited[iNodeId] = GRAY;
 
		for (auto& it : vecGraph[iNodeId])
		{
			int iNeighbourNodeId = it.first;
 
			if ((vecVisited[iNeighbourNodeId] == GRAY)
				|| (vecVisited[iNeighbourNodeId] == WHITE && cyclic_directed(cyclic_directed, iNeighbourNodeId)))
			{
				return true;
			}
		}
 
		vecVisited[iNodeId] = BLACK;
		return false;
	};
 
	auto cyclic_undirected = [&](auto&& cyclic_undirected, int iNodeId, int iParentNodeId) -> bool
	{
		vecVisited[iNodeId] = BLACK;
 
		for (auto& it : vecGraph[iNodeId])
		{
			int iNeighbourNodeId = it.first;
 
			if (vecVisited[iNeighbourNodeId] == WHITE)
			{
				if (cyclic_undirected(cyclic_undirected, iNeighbourNodeId, iNodeId))
				{
					return true;
				}
			}
			else if (iNeighbourNodeId != iParentNodeId)
			{
				return true;
			}
		}
		return false;
	};
 
	for (unsigned int i = 0; i < vecGraph.size(); ++i)
	{
		if (vecVisited[i] == WHITE)
		{
			if (m_fDirected)
			{
				if (cyclic_directed(cyclic_directed, i))
				{
					return true;
				}
			}
			else if (cyclic_undirected(cyclic_undirected, i, -1))
			{
				return true;
			}
		}
	}
 
	return false;
}
 
using namespace std;
 
int main()
{
    ios_base::sync_with_stdio(0);
    int testCases,vertices,edges,startNode,endNode;
    long long weight;
    cin>>testCases;
    while(testCases--)
    {
        cin>>vertices>>edges;
        vector<pair<int,int>> vecEdges;
        vector<long long> vecWeights;
        for(int i = 0; i < edges; ++i)
        {
            cin>>startNode>>endNode>>weight;
            --startNode;
            --endNode;
            vecEdges.push_back({startNode, endNode});
            vecWeights.push_back(weight);
        }
        cin>>startNode>>endNode;
        --startNode;
        --endNode;
        Graph<int> g = Graph<int> (vecEdges, vecWeights, vertices /**graph size*/, true /**directed graph*/);
        long long result = g.dijkstra(startNode)[endNode];
        if (result != LLONG_MAX)
        {
            cout<<result<<"\n";
        }
        else
        {
            cout<<"NO\n";
        }
    }
    return 0;
}
 
