//
// Created by yumei Hu on 2023/12/19.
//

#include "Contraction.h"

using namespace szx;

Contraction::Contraction(DFeedbackVertexSet input, AdjList reverseAdjList) {
    this->n = input.nodeNum;
    this->m = input.arcNum;
    this->cnt = 0;
    this->cntb = 0;
    this->edge = input.adjList;
    this->subtree.resize(input.nodeNum);
    vector<bool> vec(input.nodeNum, false);
    this->instack = vec;
    this->dfn.resize(input.nodeNum);
    this->low.resize(input.nodeNum);
    this->indegree.resize(input.nodeNum);
    this->outdegree.resize(input.nodeNum);
    this->reverseAdjList = reverseAdjList;
}

Contraction::~Contraction() {

}

void Contraction::domeOperation() {

}

void Contraction:: pieOperation() {
    //第一种情况，删除一个scc到另一个scc的边
}

//五步化简法
void Contraction::preContraction() {
    calculateInAndOutDegree();
    queue<int> queue;
    vector<bool> visited(this->indegree.size(), false);
    updateVisitedList(queue, visited);

    while (!queue.empty()) {
        int v = queue.front();
        queue.pop();
        if (this->indegree[v] == 0) {//IN0, 删除入度为0的顶点及其关联边
           deleteNode(v);
        } else if (this->outdegree[v] == 0) {//OUT0，删除出度为0的顶点及其关联边
            deleteNode(v);
        } else if (isSelfLoop(v)) {//LOOP
            this->cutset.push_back(v);
            deleteNode(v);
        } else if (this->indegree[v] == 1) {//IN1
            int u = this->reverseAdjList[v][0];//u是v的唯一前序，把v合并到u
            combineForward(v, u);
            deleteNode(v);
        } else if (this->outdegree[v] == 1) {//OUT1
            int u = this->edge[v][0];//u是v的唯一后序，把v合并到u
            combineAfterward(v, u);
            deleteNode(v);
        }
        updateVisitedList(queue, visited);
    }
}

bool Contraction::isSelfLoop(int nodeId){
    if (this->edge[nodeId].size() > 0) {
        vector<int>::iterator it = std::find(this->edge[nodeId].begin(), this->edge[nodeId].end(), nodeId);
        if (it != this->edge[nodeId].end()) {
            return true;
        }
    }
    return false;
}

void Contraction::updateVisitedList(queue<int>& queue, vector<bool>& visited) {
    for (int i = 0; i < this->indegree.size(); i++) {
        if (!visited[i]) {
            if (this->indegree[i] == 0 && this->outdegree[i] == 0) {
                visited[i] = true;
            } else if (this->indegree[i] <= 1 || this->outdegree[i] <= 1) {
                queue.push(i);
                visited[i] = true;
            }
        }
    }
}

void Contraction::combineForward(int v, int u) {
    vector<int> toNodes = this->edge[v];
    for (int toNode : toNodes) {
        this->edge[u].push_back(toNode);
        this->reverseAdjList[toNode].push_back(u);
        addOutdegree(u);
        addIndegree(toNode);
    }
}

void Contraction::combineAfterward(int v, int u) {
    vector<int> fromNodes = this->reverseAdjList[v];
    for (int fromNode : fromNodes) {
        this->edge[fromNode].push_back(u);
        this->reverseAdjList[u].push_back(fromNode);
        addIndegree(u);
        addOutdegree(fromNode);
    }
}

void Contraction::deleteNode(int nodeId) {
    vector<int> fromNodes = this->reverseAdjList.at(nodeId);
    vector<int> toNodes = this->edge.at(nodeId);
    for (int fromNode : fromNodes) {
        this->edge[fromNode].erase(std::remove(this->edge[fromNode].begin(), this->edge[fromNode].end(), nodeId),
                            this->edge[fromNode].end());
        minusOutdegree(fromNode);
    }
    for (int toNode : toNodes) {
        this->reverseAdjList[toNode].erase(std::remove(this->reverseAdjList[toNode].begin(),this->reverseAdjList[toNode].end(), nodeId),
                                      this->reverseAdjList[toNode].end());
        minusIndegree(toNode);
    }
    this->edge[nodeId].clear();
    this->reverseAdjList[nodeId].clear();
    clearIndegree(nodeId);
    clearOutdegree(nodeId);
}

void Contraction::calculateInAndOutDegree() {
    //计算每个顶点的入度和出度，分别保存在indegree、outdegree数组中
    for (int i = 0; i < this->n; i++) {
        this->outdegree[i] = this->edge[i].size();
        for (int j = 0; j < this->edge[i].size(); j++) {
            ++this->indegree[edge[i][j]];
        }
    }
}

void Contraction::tarjan(int u) {
    ++cnt;
    dfn[u] = low[u] = cnt;
    s.push(u);
    instack[u] = true;
    for (int i=0; i < edge[u].size(); i++) {
        int v = edge[u][i];
        if (!dfn[v]) {
            tarjan(v);
            low[u] = min(low[u],low[v]);
        }
        else if (instack[v]) {
            low[u] = min(low[u],dfn[v]);
        }
    }
    if (dfn[u] == low[u]) {
        ++cntb;
        int node;
        do {
            node=s.top();
            s.pop();
            instack[node]=false;
            subtree[cntb - 1].push_back(node);
        } while (node != u);
    }
}