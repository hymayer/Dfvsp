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
    queue<int> queue;
    vector<bool> inqueue(indegree.size(), false);
    for (int i = 0; i < this->indegree.size(); i++) {
        if (this->indegree[i] <= 1 || this->outdegree[i] <= 1) {
            queue.push(i);
            inqueue[i] = true;
        }
    }

    while (!queue.empty()) {
        int v = queue.front();
        queue.pop();
        //IN0, 删除入度为0的顶点及其关联边
        if (this->indegree[v] == 0) {
            this->dagV.push_back(v);
            for (int j : this->edge[v]) {
                --this->indegree[j];
                if (!inqueue[j] && this->indegree[j] <= 1) {
                    queue.push(j);
                    inqueue[j] = true;
                }
                this->reverseAdjList[j].erase(std::remove(this->reverseAdjList[j].begin(),this->reverseAdjList[j].end(),v),
                                              this->reverseAdjList[j].end());
            }
            this->edge[v].clear();
        }
        //OUT0，删除出度为0的顶点及其关联边
        else if (this->outdegree[v] == 0) {
            this->dagV.push_back(v);
            for (int j : this->reverseAdjList[v]) {
                --this->outdegree[j];
                if (!inqueue[j] && this->outdegree[j] <= 1) {
                    queue.push(j);
                    inqueue[j] = true;
                }
                this->edge[j].erase(std::remove(this->edge[j].begin(),this->edge[j].end(), v),
                                    this->edge[j].end());
            }
            this->reverseAdjList[v].clear();
        }
        //LOOP
        else if (std::find(this->edge[v].begin(), this->edge[v].end(), v) != this->edge[v].end()) {
            this->cutset.push_back(v);
            this->indegree[v] = 0;
            this->outdegree[v] = 0;
            for (int i : this->edge[v]) {
                --this->indegree[i];
                if (!inqueue[i] && this->indegree[i] <= 1) {
                    queue.push(i);
                    inqueue[i] = true;
                }
            }
            this->edge[v].clear();
            for (int j : this->reverseAdjList[v]) {
                --this->outdegree[j];
                if (!inqueue[j] && this->outdegree[j] <= 1) {
                    queue.push(j);
                    inqueue[j] = true;
                }
            }
            this->reverseAdjList[v].clear();
        }
        //IN1
        else if (this->indegree[v] == 1) {
            int u = this->reverseAdjList[v][0];//u是v的唯一前序，把v合并到u
            this->indegree[v] = 0;
            this->outdegree[v] = 0;
            --this->outdegree[u];
            for (int i : this->edge[v]) {
                this->edge[u].push_back(i);
                ++this->outdegree[u];
                this->reverseAdjList[i].erase(std::remove(this->reverseAdjList[i].begin(),this->reverseAdjList[i].end(), v),
                                              this->reverseAdjList[i].end());
                this->reverseAdjList[i].push_back(u);
            }
            this->reverseAdjList[v].clear();
            this->edge[v].clear();
            this->edge[u].erase(std::remove(this->edge[u].begin(), this->edge[u].end(), v), this->edge[u].end());
        }
        //OUT1
        else if (this->outdegree[v] == 1) {
            int u = this->edge[v][0];//u是v的唯一后序，把v合并到u
            this->indegree[v] = 0;
            this->outdegree[v] = 0;
            --this->indegree[u];
            for (int i : this->reverseAdjList[v]) {
                this->edge[i].push_back(u);//所有v前序顶点的出边加上u
                ++this->indegree[u];
                this->edge[i].erase(std::remove(this->edge[i].begin(), this->edge[i].end(), v),
                                    this->edge[i].end());//所有v前序顶点的出边删除v
                this->reverseAdjList[u].push_back(i);
            }
            this->reverseAdjList[v].clear();
            this->edge[v].clear();
            this->reverseAdjList[u].erase(std::remove(this->reverseAdjList[u].begin(), this->reverseAdjList[u].end(), v),
                                          this->reverseAdjList[u].end());
        }
    }
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