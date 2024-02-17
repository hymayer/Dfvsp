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

#define MAXVEX 100//��󶥵���
#define INFINITY 65535//�������� 
#define MAXEDGE 1000

using namespace std;

namespace szx {
    //����
    struct ArcNode {
        ArcNode* next; //��һ�������ı�
        int adjvex;   //���满β�����ڶ�����е��±�
    };
    struct Vnode {
        string data; //��������
        ArcNode* firstarc; //��һ�������ڸö���ı�
    };

    class Graph_DG {
    private:
        int vexnum; //ͼ�Ķ�����
        int edge;   //ͼ�ı���
        int* indegree; //ÿ�������������
        Vnode* arc; //�ڽӱ�
    public:
        Graph_DG(int, int);
        ~Graph_DG();
        //����һ��ͼ
        void createGraph(DFeedbackVertexSet& input);
        //��ӡ�ڽӱ�
        void print();
        //������������,Kahn�㷨
        bool topological_sort();
    };

}
