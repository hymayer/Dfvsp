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
        vector<int> indegree;
        vector<int> outdegree;
        vector<int> cutset;
        vector<int> dagV;
        AdjList reverseAdjList;
    public:
        Contraction(DFeedbackVertexSet input, AdjList reverseAdjList);
        ~Contraction();
        //计算图的scc，使用tarjan算法
        void tarjan(int u);
        void pieOperation();
        void domeOperation();
        void preContraction();
        //计算每个顶点的入度和出度
        void calculateInAndOutDegree();
    };
}


#endif //DFVSP_CONTRACTION_H
