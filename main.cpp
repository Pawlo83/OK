#include <iostream>
#include <ctime>
#include "nlohmann/json.hpp"
#include <fstream>

#include <vector>
#include <random>
#include <cmath>
#include <cstdlib>

#include <queue>

#include <chrono>

using namespace std;
using json = nlohmann::json;
int x =0;
// Define random number generator
random_device rd;
mt19937 gen(rd());

void dfs(const vector<vector<int>>& matrix, int vertex, vector<bool>& visited){
    visited[vertex]=true;
    for(int i=0; i<matrix.size();++i){
        if(matrix[vertex][i]!=numeric_limits<int>::max() && !visited[i]){
            dfs(matrix,i,visited);
        }
    }
}
bool is_connected(const vector<vector<int>>& matrix){
    int n=matrix.size();
    vector<bool> visited(n,false);
    dfs(matrix,0,visited);
    // Check if all vertices were visited
    bool allVerticesVisited=true;
    for(bool v : visited){
        if(!v){
            allVerticesVisited=false;
            break;
        }
    }
    return allVerticesVisited;
}

vector<vector<int>> generateGraph(int n, float d, int r){
    int rate=int(d*(n*(n-1)/2)); // Calculate graph saturation factor
    vector<vector<int>> matrix(n,vector<int>(n,numeric_limits<int>::max())); //Create the matrix filed with infinity
    for(int i=0;i<n;i++){ // Fill the diagonal with 0s
        matrix[i][i]=0;
    }
    for(int _=0;_<rate;_++){
        int v1=uniform_int_distribution<>(0,n-1)(gen);
        int v2=uniform_int_distribution<>(0,n-1)(gen);
        while(matrix[v1][v2] != numeric_limits<int>::max() || matrix[v2][v1]!=numeric_limits<int>::max() || v1==v2){
            v1=uniform_int_distribution<>(0,n-1)(gen);
            v2=uniform_int_distribution<>(0,n-1)(gen);
        }
        uniform_int_distribution<> dis(1,r); // Set range of possible distances between vertices
        int value=dis(gen);
        matrix[v1][v2]=value;
        matrix[v2][v1]=value;
    }
    if(is_connected(matrix)==1){
   		return matrix;
    }
    else{
    	cout << "TRY AGAIN" << endl;
    	return generateGraph(n,d,r);
    }
}

void printMatrix(const vector<vector<int>>& matrix){
    int n=matrix.size();
    for(int i=0;i<n;++i){
        for (int j=0;j<n;++j){
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

void sortMatrixByColumn(vector<vector<float>>& matrix, int column) {
    sort(matrix.begin(), matrix.end(), [column](const auto& a, const auto& b) {
        return a[column] > b[column];
    });
}

void generateGraphImage(const vector<vector<int>>& matrix, const string& filename){
    ofstream dotFile("graph.dot");
    dotFile << "graph G {\n";

    int n=matrix.size();

    for(int i=0;i<n;++i){
        for(int j=i+1;j<n;++j){
            if(matrix[i][j]!=numeric_limits<int>::max()){
                dotFile << "  " << i << " -- " << j << " [label=\"" << matrix[i][j] << "\"]\n";
            }
        }
    }

    dotFile << "}\n";
    dotFile.close();

    string command="dot -Tpng graph.dot -o "+filename;
    system(command.c_str());
}

// Generate a flow matrix based on a graph
vector<vector<int>> generateFlow(const vector<vector<int>>& graph, int f){
    int n=graph.size();
    vector<vector<int>> flow_matrix(n,vector<int>(n,0));

    // Assign random flow values to edges
    for (int i=0;i<n;++i){
        for (int j=i+1;j<n;++j){
            if(graph[i][j]!=numeric_limits<int>::max()){
                flow_matrix[i][j]=uniform_int_distribution<>(1,f)(gen);
                flow_matrix[j][i]=flow_matrix[i][j];
            }
        }
    }
    return flow_matrix;
}

// Create a Graphviz DOT and PNG file for a flow matrix
void generateFlowImage(const vector<vector<int>>& flow_matrix, const string& filename){
    ofstream dotFile("flow.dot");
    dotFile << "digraph G {\n";

    int n=flow_matrix.size();

    for(int i=0;i<n;++i){
        for(int j=0;j<n;++j){
            if(flow_matrix[i][j]!=0){
                dotFile << "  " << i << " -> " << j << " [label=\"" << flow_matrix[i][j] << "\"]\n";
            }
        }
    }

    dotFile << "}\n";
    dotFile.close();

    string command="dot -Tpng flow.dot -o " + filename;
    system(command.c_str());
}

bool bfs(const vector<vector<int>>& residual, vector<int>& parent, int source, int sink){
    int n=residual.size();
    vector<bool> visited(n,false);
    queue<int> q;
    q.push(source);
    visited[source]=true;
    parent[source]=-1;
    while(!q.empty()){
        int u=q.front();
        q.pop();
        for(int v=0;v<n;++v){
            if(!visited[v] && residual[u][v]>0){
                q.push(v);
                parent[v]=u;
                visited[v]=true;
                if (v==sink)
                    return true;
            }
        }
    }
    return false;
}

void generateKarpImage(const vector<vector<int>>& flow_matrix, const vector<vector<int>>& residual, const string& filename, int iteration){
    ofstream dotFile("karp.dot");
    dotFile << "digraph G {\n";

    int n=flow_matrix.size();

    for(int i=0;i<n;++i){
        for(int j=0;j<n;++j){
            if(flow_matrix[i][j]!=0){
                dotFile << "  " << i << " -> " << j << " [label=\"" << flow_matrix[i][j] << ", " << flow_matrix[i][j]-residual[i][j] << "\", ";
                // Add color for used paths
                if (residual[i][j] != flow_matrix[i][j]){
                    dotFile << "color=\"red\"";
                }
                dotFile << "]\n";
            }
        }
    }

    dotFile << "}\n";
    dotFile.close();

    string command = "dot -Tpng karp.dot -o "+filename+to_string(iteration)+".png";
    system(command.c_str());
}

int edmonds_karp(const vector<vector<int>>& flow_matrix, int source, int sink, int iteration=-1){
    int n=flow_matrix.size();
    vector<vector<int>> residual(flow_matrix);

    int max_flow=0;
    vector<int> parent(n);

    while (bfs(residual,parent,source,sink)){
        int path_flow=numeric_limits<int>::max();

        // Find the minimum capacity along the path
        for(int v=sink;v!=source;v=parent[v]){
            int u=parent[v];
            path_flow=min(path_flow,residual[u][v]);
        }

        // Update the residual graph and flow
        for(int v= sink;v!=source;v=parent[v]){
            int u=parent[v];
            residual[u][v]-=path_flow;
            residual[v][u]+=path_flow;
		}
        max_flow+=path_flow;
    }

    // Generate graphical representation of the residual flow matrix
    if(max_flow>0){
    	generateKarpImage(flow_matrix,residual,"karp", iteration);
    	//printMatrix(residual);
    }	
    
    return max_flow;
}

vector<vector<float>> finding_single_connections(const vector<vector<int>>& graph,const vector<vector<int>>& flow_matrix, int source, int sink, const string finding_method, vector<int> iff = {}, int iff1=0){
	int n=flow_matrix.size();
	vector<vector<float>> possible_edges(n,vector<float>(4,0));
	vector<vector<int>> work_matrix(flow_matrix);
	int j=0;
	
	for(int i=0;i<n;++i){
		auto ifff=find(iff.begin(), iff.end(), i);
		if(ifff==iff.end()){
			work_matrix[sink][i]=0;
			work_matrix[i][sink]=0;
		}
	}
	//printMatrix(work_matrix);
	//cout << "^ fsc1" << endl;
	for(int i=0;i<n;++i){
		work_matrix[sink][i]=flow_matrix[sink][i];
		work_matrix[i][sink]=flow_matrix[i][sink];
		//printMatrix(work_matrix);
		//cout << "^ fsc" << endl;
		int flow;
		if(finding_method=="edmonds_karp"){
           	flow=edmonds_karp(work_matrix,source,sink,x);
           	x++;
	    }
	    //else if(finding_method=="pushRelabelMaxFlow"){
	     	//flow=pushRelabelMaxFlow(work_matrix,source,sink);
	    //}
	    else{
	     	cerr << "Invalid flow finding method." << endl;
	     	return {};
     	}
     	//cout << flow << " " << iff1 <<endl;
		if(flow-iff1>0){
			possible_edges[j][0]=i;
			possible_edges[j][1]=flow-iff1;
			possible_edges[j][2]=graph[sink][i];
			possible_edges[j][3]=static_cast<float>(flow)/graph[sink][i];
			j++;
			cout << "Max flow: " << flow-iff1 << ", edge: " << i << "," << sink << ", edge length: " << graph[sink][i] <<endl;
		}

		auto ifff=find(iff.begin(), iff.end(), i);
		if(ifff==iff.end()){
			work_matrix[sink][i]=0;
			work_matrix[i][sink]=0;
		}
		
	}
	sortMatrixByColumn(possible_edges,3);
	return possible_edges;
}

void choose_edges(const vector<vector<int>>& graph, const vector<vector<int>>& flow_matrix, vector<vector<float>> possible_edges, int source, int sink, int max_flow, const string finding_method){
	int needed_flow_inp;
	int n=possible_edges.size();
	vector<vector<int>> work_matrix(flow_matrix);
	for(int i=0;i<n;++i){
		work_matrix[sink][i]=0;
		work_matrix[i][sink]=0;
	}
	vector<int> used_edges(n,-1);
	cout << "Input needed flow: ";
	cin >> needed_flow_inp;
	int needed_flow = needed_flow_inp;

	if(needed_flow>max_flow){
		cout << "Max flow = " << max_flow << ", so its imposible to get flow that = " << needed_flow_inp << endl;
		return;
	}
	int flowFull=0;
	int i=0;
	//for(int i=0;i<n;++i){
	while(flowFull<needed_flow){
		//printMatrix(work_matrix);
		//cout << endl;
		work_matrix[sink][possible_edges[0][0]]=flow_matrix[sink][possible_edges[0][0]];
		work_matrix[possible_edges[0][0]][sink]=flow_matrix[possible_edges[0][0]][sink];
		used_edges[i]=possible_edges[0][0];
	    //for(int z=0;z<=i;z++){
	    //	cout << "cc" <<endl;
	    //	cout << used_edges[z] <<"cc"<<endl;
	    //}
		//printMatrix(work_matrix);
		//cout << endl;
		if(finding_method=="edmonds_karp"){
           	flowFull=edmonds_karp(work_matrix,source,sink,-i-2);
           	cout << flowFull << endl;
	    }

		//if(flowFull>=needed_flow){
		//	break;
		//}
		possible_edges=finding_single_connections(graph,flow_matrix,source,sink,"edmonds_karp",used_edges,flowFull);
		i++;
	}

	int full_added_lenght=0;
	cout << "Edges needed to be used:" << endl;
	for(int i=0;i<n;i++){
		if(used_edges[i]==-1){
			break;
		}
		if(i!=0){
			cout << "; ";
		}
		cout << used_edges[i] << "," << sink;
		full_added_lenght+=graph[used_edges[i]][sink];
	}
	cout<<endl;
	cout << "Full added length: "<<full_added_lenght<<endl;
}

int main(){
    int n=5; // Set the number of vertices
    int r=1000; // Set max distance
    int f=30; // Set max flow
    float d=0.9; // Set saturation
    vector<vector<int>> graph=generateGraph(n,d,r);
    vector<vector<int>> flow_matrix=generateFlow(graph,f);

    // Print generated graphs and view them as images
    //printMatrix(graph);
    //printMatrix(flow_matrix);
    generateGraphImage(graph, "graph.png");
    generateFlowImage(flow_matrix, "flow.png");

	// Set random source and sink
    int source=uniform_int_distribution<>(0, n-1)(gen);
    int sink=source;
    while(sink==source){
    	sink=uniform_int_distribution<>(0, n-1)(gen);
    }
	cout << "Source: " << source << " " << "Sink:" << sink << endl;
	
 	auto start_timeEK = chrono::steady_clock::now();
	int max_flowEK=edmonds_karp(flow_matrix,source,sink);
	cout << "Max flow: " << max_flowEK << endl;
	vector<vector<float>> possible_edgesEK=finding_single_connections(graph,flow_matrix,source,sink,"edmonds_karp");
	auto end_timeEK = chrono::steady_clock::now(); 
	auto timeEK=chrono::duration_cast<chrono::milliseconds>(end_timeEK - start_timeEK);
	cout << "Time for EK: " << timeEK.count() << endl;
	
	choose_edges(graph, flow_matrix, possible_edgesEK, source, sink, max_flowEK, "edmonds_karp");
	
    return 0;
}
