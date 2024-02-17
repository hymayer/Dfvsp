#pragma once
#include <iostream>
#include <sstream>
#include<stdio.h>
#include<queue>
#include<stack>
#include<stdlib.h>
#include<string>
#include<stack>
#include "DFeedbackVertexSet.h"

#define MAXVEX 100//最大顶点数
#define INFINITY 65535//代表无穷 
#define MAXEDGE 1000

using namespace std;

namespace szx {
    //表结点
    struct ArcNode {
        ArcNode* next; //下一个关联的边
        int adjvex;   //保存弧尾顶点在顶点表中的下标
    };
    struct Vnode {
        string data; //顶点名称
        ArcNode* firstarc; //第一个依附在该顶点的边
    };

    class Graph_DG {
    private:
        int vexnum; //图的顶点数
        int edge;   //图的边数
        int* indegree; //每个顶点的入度情况
        Vnode* arc; //邻接表
    public:
        Graph_DG(int, int);
        ~Graph_DG();
        //创建一个图
        void createGraph(DFeedbackVertexSet& input);
        //打印邻接表
        void print();
        //进行拓扑排序,Kahn算法
        bool topological_sort();
    };

}
