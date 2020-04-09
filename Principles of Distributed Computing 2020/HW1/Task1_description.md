**1. Spanning Trees**

This question considers the construction of a spanning tree in a distributed setting, where nodes are not allowed to send messages. Every node can only maintain a binary value for each incident edge, which is initially set to 0. Instead of communicating, each node can observe the XOR result of the two values for each incident edge. For example, for an edge (u,v), if node u assigns 1 and node v also assigns 1, then both nodes u and v can observe the XOR value 1⊕1=0 on the edge (u,v) . 

Given an undirected graph G=(V,E) on n nodes, the task is to develop a deterministic algorithm for constructing a spanning tree in this message-free setting. Each node has a unique ID 0,1,…,n−1. We assume that each node knows its own ID, but not the ID of its neighbors. Moreover, each node can differentiate between its incident edges, e.g., each node can label its edges locally. At the end of the constructing, the edges in the spanning tree output 1 and others edges output 0.

In the first part of the question, we consider the synchronous model: In each round, every node can observe the XOR result of each incident edge, perform some local computation, and potentially update its binary value for all of its edges. Any algorithm with a running time of O(n) will be given full points.

1. Design a synchronous algorithm that can construct a spanning tree in G. [5]
2. Design a synchronous algorithm that can construct a spanning tree in G and, in addition, one node knows when the construction is finished. [5]

Next, we consider a sequential round-robin model: Nodes are ordered according to some predefined sequence. In each round, the next node is chosen according to this sequence, and it can execute the same steps as in the synchronous setting. Once we reach the end of the sequence, we continue with the first node of the sequence.

3. Design an algorithm that constructs a spanning tree in G in the sequential round-robin model.  Analyze the round complexity of your algorithm in the best possible graph under the best-case ordering, as well as in the worst possible graph under the worst-case ordering. [10]

**2. k-max-Coloring on Trees**

In this question, we consider deterministic algorithms for coloring trees in a message-free setting. Given a tree of n nodes, we want each node to choose a color (represented by a positive integer) different from the colors of its neighbors. We say that an algorithm provides a k-max-coloring if no node has a color larger than k when the algorithm terminates. In our model, nodes operate in synchronous rounds. Initially, all nodes have color 1. In each round, each node can execute the following steps:

Observe the value of the colors of its neighbors.
Do some local computations
Possibly increase its color to an arbitrary larger integer. Note that nodes are not allowed to decrease their color.
We assume that each node has a unique ID 0,1,…,n−1 and that each node knows n and its own ID, but not the ID of its neighbors.

Under the model assumptions described above, answer the following questions:
Give an algorithm providing a n-max-coloring in O(1) rounds. [5]
Is there an algorithm providing a 3-max-coloring in O(n) rounds? [5]

**3. Maximal Trianglefree Subgraph**

Consider a network graph G=(V,E) with n vertices, each representing a computer, with unique identifiers in 1,...,n, where per round each computer can send a message to each of its neighbors. Suppose that Δ denotes the maximum degree in G and the values n and Δ are known to all computers. 


The objective is to compute a set S⊆V of vertices such that (1) the subgraph G[S] induced by S --- that is, the subgraph made of vertices in S and edges of G for which both endpoints are in S --- does not contain a triangle, i.e., three vertices that are all connected to each other, and (2) we cannot add anymore vertex from V∖S to S without violating condition (1).  


Devise a deterministic distributed algorithm that computes such a set S in f(Δ)+O(log∗n) rounds, for the smallest function f(Δ) that you can find. Each computer should know whether it's node is in S or not.

**4. Marking vertices**
Suppose somebody gives you a deterministic distributed graph algorithm that for any given graph G marks some subset S of its vertices in such a way that no two marked vertices are connected by an edge and for any vertex v∈G there is a marked vertex s∈S such that the distance of v and s is at most 100. 

Prove that there is a network graph with n nodes where the round complexity of such an algorithm must be at least Ω(log∗n).
Is there an algorithm providing a 2-max-coloring in O(n) rounds? [5]
Give an algorithm that allows each node to communicate its ID to all its neighbors in O(logn) rounds such that no node has a color larger than ⌈logn⌉+1 at the end of the algorithm. Note that we do not require nodes to have a different color from their neighbors. [5]
Give an algorithm providing a O(logcn)-max-coloring in O(logc+1n) rounds for some c∈N. [10]
