//
// Created by yumei Hu on 2023/12/19.
//

#ifndef DFVSP_CONTRACTION_H
#define DFVSP_CONTRACTION_H

#include <vector>
#include <stack>
#include <queue>
#include "DFeedbackVertexSet.h"

using namespace std;

namespace szx {
    class Contraction {
    private:
        int n, m, cnt, cntb;
        vector<vector<int>> edge;
        vector<vector<int>> subtree;
        vector<bool> instack;
        vector<int> dfn, low;
        stack<int> s;
        vector<int> indegree; //记录所有顶点的入度
        vector<int> outdegree; //记录所有顶点的出度
        vector<int> cutset;
        AdjList reverseAdjList;
    public:
        Contraction(DFeedbackVertexSet& input, AdjList& reverseAdjList);
        ~Contraction();
        //计算图的scc，使用tarjan算法
        void tarjan(int u);
        void pieOperation();
        void domeOperation();
        void LLContraction();
        //计算每个顶点的入度和出度
        void calculateInAndOutDegree();
        //删除所有点i关联的边
        void deleteNode(int i);
        //合并顶点v到前序顶点u
        void combineForward(int v, int u);
        //合并顶点v到后序顶点u
        void combineAfterward(int v, int u);
        //更新队列
        void updateVisitedList(queue<int>& queue, vector<bool>& visited);
        //判断当前节点是否自循环
        bool isSelfLoop(int nodeId);
        vector<int> getCutset() { return this->cutset; };
        int getIndegree(int nodeId) { return this->indegree[nodeId]; }
        void addIndegree(int nodeId) { this->indegree[nodeId]++; }
        void minusIndegree(int nodeId) { this->indegree[nodeId]--; }
        void clearIndegree(int nodeId) { this->indegree[nodeId] = 0; }
        int getOutdegree(int nodeId) { return  this->outdegree[nodeId]; }
        void addOutdegree(int nodeId) { this->outdegree[nodeId]++; }
        void minusOutdegree(int nodeId) { this->outdegree[nodeId]--; }
        void clearOutdegree(int nodeId) { this->outdegree[nodeId] = 0; }
    };
}


#endif //DFVSP_CONTRACTION_H
